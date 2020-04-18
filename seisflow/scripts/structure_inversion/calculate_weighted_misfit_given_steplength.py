"""
calculate_weighted_misfit_given_steplength: calculate the misfit given a step length perturbation.
"""
from glob import glob
from os.path import basename, join

import numpy as np
from mpi4py import MPI

from ...tasks.line_search.line_search_structure import calculate_weighted_misfit, get_perturbed_vir_sync
from ...utils.asdf_io import VirAsdf
from ...utils.load_files import load_first_arrival_baz_evdp, load_windows
from ...setting import SURFACE_THRESHOLD, LINE_SEARCH_PERTURBATION
from ...utils.get_path import get_data_asdf_fnames, get_sync_asdf_fnames

comm = MPI.COMM_WORLD  # pylint: disable=c-extension-no-member
size = comm.Get_size()
rank = comm.Get_rank()


def get_event_this_rank(cmts_directory):
    """
    get_event_this_rank: get the event used in this rank.
    """
    all_cmts_fpath = sorted(glob(join(cmts_directory, "*")))
    assert len(all_cmts_fpath) == size
    cmt_fpath_this_rank = all_cmts_fpath[rank]
    cmt_this_rank = basename(cmt_fpath_this_rank)
    return cmt_this_rank


def get_consider_surface(cmt_this_rank, evdp):
    """
    Get the status that if we should consider the surface wave.
    """
    rep_net_sta = list(evdp[cmt_this_rank].keys())[0]
    if (evdp[cmt_this_rank][rep_net_sta] <= SURFACE_THRESHOLD):
        return True
    else:
        return False


def load_virasdf_files(
        cmt_this_rank, min_periods, max_periods, data_asdf_directory, sync_raw_directory, sync_perturbed_directory):
    """
    load_virasdf_files: get the data and sync virasdf files for the surface/body waves.
    """
    # data_asdf_body_path, sync_asdf_body_path, data_asdf_surface_path, sync_asdf_surface_path = get_asdf_fnames(
    #     cmt_this_rank, min_periods, max_periods, data_asdf_directory, sync_asdf_directory)
    data_asdf_body_path, data_asdf_surface_path = get_data_asdf_fnames(
        cmt_this_rank, min_periods, max_periods, data_asdf_directory)
    sync_asdf_body_path_raw, sync_asdf_surface_path_raw = get_sync_asdf_fnames(
        cmt_this_rank, min_periods, max_periods, sync_raw_directory)
    sync_asdf_body_path_perturbed, sync_asdf_surface_path_perturbed = get_sync_asdf_fnames(
        cmt_this_rank, min_periods, max_periods, sync_perturbed_directory)

    # * load virasdf files
    data_virasdf_body = VirAsdf()
    sync_virasdf_body_raw = VirAsdf()
    sync_virasdf_body_perturbed = VirAsdf()
    data_virasdf_surface = VirAsdf()
    sync_virasdf_surface_raw = VirAsdf()
    sync_virasdf_surface_perturbed = VirAsdf()
    data_virasdf_body.read_asdf(data_asdf_body_path)
    sync_virasdf_body_raw.read_asdf(sync_asdf_body_path_raw)
    sync_virasdf_body_perturbed.read_asdf(sync_asdf_body_path_perturbed)
    data_virasdf_surface.read_asdf(data_asdf_surface_path)
    sync_virasdf_surface_raw.read_asdf(sync_asdf_surface_path_raw)
    sync_virasdf_surface_perturbed.read_asdf(sync_asdf_surface_path_perturbed)

    return data_virasdf_body, sync_virasdf_body_raw, sync_virasdf_body_perturbed, data_virasdf_surface, sync_virasdf_surface_raw, sync_virasdf_surface_perturbed


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--cmts_directory', required=True, type=str, help="the reference cmt solution directory")
    @click.option('--windows_directory', required=True, type=str, help="the windows directory")
    @click.option('--data_info_directory', required=True, type=str, help="the data info directory")
    @click.option('--data_asdf_directory', required=True, type=str, help="the processed data directory")
    @click.option('--sync_raw_directory', required=True, type=str, help="the raw processed sync directory")
    @click.option('--sync_perturbed_directory', required=True, type=str, help="the perturbed processed sync directory")
    @click.option('--stations_path', required=True, type=str, help="the stations path in specfem format")
    @click.option('--step_length', required=True, type=float, help=f"the step length used, assume the perturbation is {LINE_SEARCH_PERTURBATION:.3f}")
    @click.option('--min_periods', required=True, type=str, help="the min periods, i.e. min_body,min_surface")
    @click.option('--max_periods', required=True, type=str, help="the max periods, i.e. max_body,max_surface")
    def main(cmts_directory, windows_directory, data_info_directory, data_asdf_directory, sync_raw_directory, sync_perturbed_directory,
             stations_path, step_length, min_periods, max_periods):
        """
        main: calculate the weighted misfit for all the events.
        """
        # * get gcmtid
        cmt_this_rank = get_event_this_rank(cmts_directory)
        print(f"rank{rank}: cmt_this_rank{cmt_this_rank}")
        # * get windows
        windows = load_windows(cmt_this_rank, windows_directory)
        # * get first_arrival_zr, first_arrival_t, baz, evdp
        first_arrival_zr, first_arrival_t, baz, evdp = load_first_arrival_baz_evdp(
            data_info_directory)
        # * get consider_surface
        consider_surface = get_consider_surface(cmt_this_rank, evdp)
        # * get virasdf
        data_virasdf_body, sync_virasdf_body_raw, sync_virasdf_body_perturbed, data_virasdf_surface, sync_virasdf_surface_raw, sync_virasdf_surface_perturbed = load_virasdf_files(
            cmt_this_rank, min_periods, max_periods, data_asdf_directory, sync_raw_directory, sync_perturbed_directory)
        # * get perturbed virasdf for sync
        sync_virasdf_body = get_perturbed_vir_sync(
            sync_virasdf_body_raw, sync_virasdf_body_perturbed, step_length)
        sync_virasdf_surface = get_perturbed_vir_sync(
            sync_virasdf_surface_raw, sync_virasdf_surface_perturbed, step_length)
        # * get stations
        stations = np.loadtxt(stations_path, dtype=np.str)
        # * calculate the misfit
        # weighted_misfit_raw = calculate_weighted_misfit(windows, consider_surface, data_virasdf_body, sync_virasdf_body_raw, data_virasdf_surface, sync_virasdf_surface_raw, first_arrival_zr, first_arrival_t, baz,
        #                                                 stations)
        weighted_misfit = calculate_weighted_misfit(windows, consider_surface, data_virasdf_body, sync_virasdf_body, data_virasdf_surface, sync_virasdf_surface, first_arrival_zr, first_arrival_t, baz,
                                                    stations)
        # * print the result
        if(rank == 0):
            print("=" * 20)
            print(f"step length: {step_length}")
            # print(
            #     f"weighted misfit for all events[raw]: {weighted_misfit_raw}")
            print(f"weighted misfit for all events: {weighted_misfit}")
            print("=" * 20)

    main()  # pylint: disable=no-value-for-parameter
