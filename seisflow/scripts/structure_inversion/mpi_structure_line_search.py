"""
mpi_structure_line_search.py: perform the structure line search in parallel for all the events.
"""
import sys
from glob import glob
from os.path import basename, join

import numpy as np
from mpi4py import MPI

from ...tasks.line_search.line_search_structure import (
    calculate_weighted_misfit, get_perturbed_vir_sync)
from ...utils.asdf_io import VirAsdf
from ...utils.get_path import get_data_asdf_fnames, get_sync_asdf_fnames
from ...utils.load_files import load_first_arrival_baz_evdp, load_windows
from ...setting import SURFACE_THRESHOLD

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


def do_search(each_search_step, windows, consider_surface, data_virasdf_body, data_virasdf_surface, sync_virasdf_body_raw, sync_virasdf_body_perturbed,
              sync_virasdf_surface_raw, sync_virasdf_surface_perturbed, stations, first_arrival_zr, first_arrival_t, baz):
    """
    do_search: calculate the perturbed waveform and do the line search.
    """
    # # * get perturbed virasdf for sync
    sync_virasdf_body = get_perturbed_vir_sync(
        sync_virasdf_body_raw, sync_virasdf_body_perturbed, each_search_step)
    sync_virasdf_surface = get_perturbed_vir_sync(
        sync_virasdf_surface_raw, sync_virasdf_surface_perturbed, each_search_step)
    # # * calculate the misfit
    weighted_misfit = calculate_weighted_misfit(windows, consider_surface, data_virasdf_body, sync_virasdf_body,
                                                data_virasdf_surface, sync_virasdf_surface, first_arrival_zr, first_arrival_t, baz, stations)
    return weighted_misfit


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
    @click.option('--min_periods', required=True, type=str, help="the min periods, i.e. min_body,min_surface")
    @click.option('--max_periods', required=True, type=str, help="the max periods, i.e. max_body,max_surface")
    @click.option('--search_range', required=True, type=str, help="the search range, eg: 0,0.03, should always from 0.")
    @click.option('--search_step', required=True, type=float, help="the searching step length, eg: 0.003, if minus searching range, -0.003")
    def main(cmts_directory, windows_directory, data_info_directory, data_asdf_directory, sync_raw_directory, sync_perturbed_directory,
             stations_path,  min_periods, max_periods, search_range, search_step):
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
        # * for each value in the searching list, do the line search
        search_range = tuple(map(float, search_range.split(",")))
        search_step_length_list = np.arange(
            search_range[0], search_range[1] + search_step, search_step)
        search_misfit_result = []
        stations = np.loadtxt(stations_path, dtype=np.str)
        for each_search_step in search_step_length_list:
            search_misfit_result.append(
                do_search(each_search_step, windows, consider_surface, data_virasdf_body, data_virasdf_surface, sync_virasdf_body_raw, sync_virasdf_body_perturbed,
                          sync_virasdf_surface_raw, sync_virasdf_surface_perturbed, stations, first_arrival_zr, first_arrival_t, baz))
            print(
                f"step: {each_search_step}; misfit: {search_misfit_result[-1]}")
            sys.stdout.flush()
        # * print the result
        if (rank == 0):
            search_misfit_result = np.array(search_misfit_result)
            print("=" * 20)
            for step_length, misfit_value in zip(search_step_length_list, search_misfit_result):
                print(f"[step length {step_length}] {misfit_value}")
            print("=" * 20)
            # * for the later model update, we should store a txt file about the step length
            optimized_step_length = search_step_length_list[np.argmin(
                search_misfit_result)]
            np.savetxt("./STEP_LENGTH", [optimized_step_length], fmt="%.3f")

    main()  # pylint: disable=no-value-for-parameter
