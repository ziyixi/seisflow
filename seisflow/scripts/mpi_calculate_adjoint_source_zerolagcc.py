"""
mpi_calculate_adjoint_source_zerolagcc.py: calculate adjoint source for multiple events in parallel using MPI.
"""
from glob import glob
from os.path import basename, join

import numpy as np
import pyasdf
from mpi4py import MPI

from ..tasks.adjoint.calculate_adjoint_source_zerolagcc_one_event import \
    calculate_adjoint_source_zerolagcc_one_event
from ..utils.asdf_io import VirAsdf
from ..utils.get_path import get_asdf_fnames
from ..utils.load_files import load_first_arrival_baz_evdp, load_pickle
from ..utils.setting import (CC_THRESHOLD, DELTAT_THRESHOLD, SNR_THRESHOLD,
                             SURFACE_THRESHOLD)

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
# * test is passed for this script on 01/07/2020


def get_gcmtid_this_rank(misfit_windows_directory):
    all_files = sorted(glob(join(misfit_windows_directory, "*pkl")))
    files_used_this_rank = np.array_split(all_files, size)[rank]
    return files_used_this_rank


def load_misfit_windows_this_rank(files_used_this_rank):
    misfit_windows_this_rank = {}
    for each_file in files_used_this_rank:
        gcmtid = basename(each_file).split(".")[0]
        misfit_windows_this_rank[gcmtid] = load_pickle(each_file)
    return misfit_windows_this_rank


def load_stations(stations_path):
    return np.loadtxt(stations_path, dtype=np.str)


def get_consider_surface(gcmtid, evdp):
    rep_net_sta = list(evdp[gcmtid].keys())[0]
    evdp_used = evdp[gcmtid][rep_net_sta]
    # evdp stored in the pickle file is already in the unit of km
    if (evdp_used <= SURFACE_THRESHOLD):
        return True
    else:
        return False


def save_adjoint_to_asdf(adjoint_source_zerolagcc, output_directory, gcmtid):
    """
    save_adjoint_to_asdf: save the adjoint source to asdf format to be read by Specfem.
    """
    output_fname = join(output_directory, f"{gcmtid}.h5")
    components = ["MXE", "MXN", "MXZ"]
    with pyasdf.ASDFDataSet(output_fname, mode="w", mpi=False) as output_asdf:
        for net_sta in adjoint_source_zerolagcc:
            adjoint_source_length = adjoint_source_zerolagcc[net_sta].shape[1]
            for index_component in range(3):
                component = components[index_component]
                specfem_adj_source = adjoint_source_zerolagcc[net_sta][index_component, :]
                tag = net_sta.replace(".", "_") + "_" + \
                    components[index_component]
                output_asdf.add_auxiliary_data(
                    data=specfem_adj_source, data_type="AdjointSources", path=tag, parameters={})


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--misfit_windows_directory', required=True, type=str, help="misfit windows directory")
    @click.option('--stations_path', required=True, type=str, help="stations path in Specfem format")
    @click.option('--raw_sync_directory', required=True, type=str, help="raw sync directory")
    @click.option('--sync_directory', required=True, type=str, help="processed sync directory")
    @click.option('--data_directory', required=True, type=str, help="processed data directory")
    @click.option('--data_info_directory', required=True, type=str, help="data info directory")
    @click.option('--output_directory', required=True, type=str, help="adjoint source asdf output directory")
    @click.option('--body_band', required=True, type=str, help="body wave filtering band, time_min,time_max")
    @click.option('--surface_band', required=True, type=str, help="surface wave filtering band, time_min,time_max")
    def main(misfit_windows_directory, stations_path, raw_sync_directory, sync_directory, data_directory,
             data_info_directory, output_directory, body_band, surface_band):
        # * prepare parameters
        # load misfit windows
        files_used_this_rank = get_gcmtid_this_rank(misfit_windows_directory)
        misfit_windows_this_rank = load_misfit_windows_this_rank(
            files_used_this_rank)
        stations = load_stations(stations_path)
        # load setting
        snr = SNR_THRESHOLD
        cc = CC_THRESHOLD
        deltat = DELTAT_THRESHOLD
        # get snr_threshold,cc_threshold,deltat_threshold,body_band,surface_band
        snr_threshold = tuple(map(float, snr.split(",")))
        cc_threshold = tuple(map(float, cc.split(",")))
        deltat_threshold = tuple(map(float, deltat.split(",")))
        body_band = tuple(map(float, body_band.split(",")))
        surface_band = tuple(map(float, surface_band.split(",")))
        # get raw_sync_asdf_path, sync_asdf_body_path, data_asdf_body_path
        # sync_asdf_surface_path, data_asdf_surface_path
        gcmtids_this_rank = [basename(item).split(".")[0]
                             for item in files_used_this_rank]
        raw_sync_asdf_paths_this_rank = [
            join(raw_sync_directory, f"{item}.h5") for item in gcmtids_this_rank]
        data_asdf_body_path_this_rank, sync_asdf_body_path_this_rank, data_asdf_surface_path_this_rank, sync_asdf_surface_path_this_rank = [], [], [], []
        for each_gcmtid in gcmtids_this_rank:
            min_periods = f"{body_band[0]},{surface_band[0]}"
            max_periods = f"{body_band[1]},{surface_band[1]}"
            data_asdf_body_path, sync_asdf_body_path, data_asdf_surface_path, sync_asdf_surface_path = get_asdf_fnames(
                each_gcmtid, min_periods, max_periods, data_directory, sync_directory)
            data_asdf_body_path_this_rank.append(data_asdf_body_path)
            sync_asdf_body_path_this_rank.append(sync_asdf_body_path)
            data_asdf_surface_path_this_rank.append(data_asdf_surface_path)
            sync_asdf_surface_path_this_rank.append(sync_asdf_surface_path)

        # * calculate the adjoint source and save to asdf format
        files_number_this_rank = len(files_used_this_rank)
        _, _, _, evdp_dict = load_first_arrival_baz_evdp(
            data_info_directory)
        for index in range(files_number_this_rank):
            gcmtid = gcmtids_this_rank[index]
            misfit_windows = misfit_windows_this_rank[gcmtid]
            raw_sync_asdf_path = raw_sync_asdf_paths_this_rank[index]
            # ! fix a serious bug here!!!
            sync_asdf_body_path = sync_asdf_body_path_this_rank[index]
            data_asdf_body_path = data_asdf_body_path_this_rank[index]
            data_asdf_surface_path = data_asdf_surface_path_this_rank[index]
            sync_asdf_surface_path = sync_asdf_surface_path_this_rank[index]
            consider_surface = get_consider_surface(gcmtid, evdp_dict)
            # build vir asdf
            raw_sync_virasdf = VirAsdf()
            raw_sync_virasdf.read_asdf(raw_sync_asdf_path)
            data_virasdf_body = VirAsdf()
            data_virasdf_body.read_asdf(data_asdf_body_path)
            sync_virasdf_body = VirAsdf()
            sync_virasdf_body.read_asdf(sync_asdf_body_path)
            data_virasdf_surface = VirAsdf()
            data_virasdf_surface.read_asdf(data_asdf_surface_path)
            sync_virasdf_surface = VirAsdf()
            sync_virasdf_surface.read_asdf(sync_asdf_surface_path)

            # run
            adjoint_source_zerolagcc = calculate_adjoint_source_zerolagcc_one_event(misfit_windows, stations, raw_sync_virasdf, snr_threshold, cc_threshold, deltat_threshold, body_band, surface_band,
                                                                                    consider_surface, sync_virasdf_body, data_virasdf_body, sync_virasdf_surface, data_virasdf_surface)
            save_adjoint_to_asdf(adjoint_source_zerolagcc,
                                 output_directory, gcmtid)
    main()
