"""
slurm_calculate_misfit_windows.py: calculate misfit windows
"""
from functools import partial
from glob import glob
from os.path import basename, join

import numpy as np
from mpi4py import MPI

from ..tasks.windows.calculate_misfit_windows import calculate_misfit_windows
from ..utils.load_files import (load_first_arrival_baz_evdp, load_pickle,
                                load_windows)
from ..utils.save_files import save_pickle_event

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


def get_asdf_fnames(gcmtid, min_periods, max_periods, data_asdf_directory, sync_asdf_directory):
    """
    get_asdf_fnames: get asdf fnames, min_periods=min_body,min_surface and the same with surface.
    """
    # we have to make sure the asdf files in asdf_directory is complete.
    body_min_period, surface_min_period = min_periods.split(",")
    body_max_period, surface_max_period = max_periods.split(",")
    data_asdf_body_path = join(
        data_asdf_directory, f"{gcmtid}.preprocessed_{body_min_period}s_to_{body_max_period}s")
    data_asdf_surface_path = join(
        data_asdf_directory, f"{gcmtid}.preprocessed_{surface_min_period}s_to_{surface_max_period}s")
    sync_asdf_body_path = join(
        sync_asdf_directory, f"{gcmtid}.preprocessed_{body_min_period}s_to_{body_max_period}s")
    sync_asdf_surface_path = join(
        sync_asdf_directory, f"{gcmtid}.preprocessed_{surface_min_period}s_to_{surface_max_period}s")

    return data_asdf_body_path, sync_asdf_body_path, data_asdf_surface_path, sync_asdf_surface_path


def kernel(gcmtid, windows_directory, min_periods, max_periods, data_asdf_directory, sync_asdf_directory, data_info_directory):
    """
    kernel to calculate misfit windows for each gcmtid (event)
    """
    windows = load_windows(gcmtid, windows_directory)
    data_asdf_body_path, sync_asdf_body_path, data_asdf_surface_path, sync_asdf_surface_path = get_asdf_fnames(
        gcmtid, min_periods, max_periods, data_asdf_directory, sync_asdf_directory)
    first_arrival_zr, first_arrival_t, baz, evdp = load_first_arrival_baz_evdp(
        data_info_directory)
    consider_surface = get_consider_surface(gcmtid, evdp)
    misfit_windows = calculate_misfit_windows(windows, consider_surface, data_asdf_body_path,
                                              sync_asdf_body_path, data_asdf_surface_path, sync_asdf_surface_path, first_arrival_zr, first_arrival_t, baz)
    return misfit_windows


def get_consider_surface(gcmtid, evdp):
    rep_net_sta = list(evdp[gcmtid].keys())[0]
    evdp_used = evdp[gcmtid][rep_net_sta]
    # evdp stored in the pickle file is already in the unit of km
    if (evdp_used <= 150):
        return True
    else:
        return False


def get_gcmtids_this_rank(windows_directory):
    all_files = glob(join(windows_directory, "*pkl"))
    all_gcmtids = [basename(item).split(".")[0] for item in all_files]
    return np.array_split(all_gcmtids, size)[rank]


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--windows_directory', required=True, type=str, help="windows directory")
    @click.option('--output_directory', required=True, type=str, help="misfit window directory to output")
    @click.option('--min_periods', required=True, type=str, help="body_wave_period,surface_wave_period")
    @click.option('--max_periods', required=True, type=str, help="body_wave_period,surface_wave_period")
    @click.option('--data_asdf_directory', required=True, type=str, help="data asdf directory")
    @click.option('--sync_asdf_directory', required=True, type=str, help="sync asdf directory")
    @click.option('--data_info_directory', required=True, type=str, help="data info directory")
    def main(windows_directory, output_directory, min_periods, max_periods, data_asdf_directory, sync_asdf_directory, data_info_directory):
        """
        Calculate misfit windows for all the gcmtids (events) in windows_directory
        """
        kernel_used = partial(
            windows_directory=windows_directory, min_periods=min_periods, max_periods=max_periods, data_asdf_directory=data_asdf_directory,
            sync_asdf_directory=sync_asdf_directory, data_info_directory=data_info_directory)
        gcmtids_this_rank = get_gcmtids_this_rank(windows_directory)
        for each_gcmtid in gcmtids_this_rank:
            misfit_windows = kernel_used(each_gcmtid)
            save_pickle_event(misfit_windows, output_directory, each_gcmtid)

    main()
