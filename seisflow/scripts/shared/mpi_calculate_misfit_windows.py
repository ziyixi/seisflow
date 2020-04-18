"""
mpi_calculate_misfit_windows.py: calculate misfit windows
"""
from functools import partial
from glob import glob
from os.path import basename, join

import numpy as np
from mpi4py import MPI

from ...tasks.windows.calculate_misfit_windows import calculate_misfit_windows
from ...utils.asdf_io import VirAsdf
from ...utils.get_path import get_asdf_fnames
from ...utils.load_files import (load_first_arrival_baz_evdp, load_windows)
from ...utils.save_files import save_pickle_event
from ...setting import SURFACE_THRESHOLD

comm = MPI.COMM_WORLD  # pylint: disable=c-extension-no-member
size = comm.Get_size()
rank = comm.Get_rank()
# * test is passed for this script on 01/07/2020


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
    # generate vir asdf
    data_virasdf_body = VirAsdf()
    data_virasdf_body.read_asdf(data_asdf_body_path)
    sync_virasdf_body = VirAsdf()
    sync_virasdf_body.read_asdf(sync_asdf_body_path)
    data_virasdf_surface = VirAsdf()
    data_virasdf_surface.read_asdf(data_asdf_surface_path)
    sync_virasdf_surface = VirAsdf()
    sync_virasdf_surface.read_asdf(sync_asdf_surface_path)

    # run
    misfit_windows = calculate_misfit_windows(windows, consider_surface, data_virasdf_body,
                                              sync_virasdf_body, data_virasdf_surface, sync_virasdf_surface, first_arrival_zr, first_arrival_t, baz)
    return misfit_windows


def get_consider_surface(gcmtid, evdp):
    rep_net_sta = list(evdp[gcmtid].keys())[0]
    evdp_used = evdp[gcmtid][rep_net_sta]
    # evdp stored in the pickle file is already in the unit of km
    if (evdp_used <= SURFACE_THRESHOLD):
        return True
    else:
        return False


def get_gcmtids_this_rank(windows_directory):
    all_files = sorted(glob(join(windows_directory, "*pkl")))
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
        kernel_used = partial(kernel,
                              windows_directory=windows_directory, min_periods=min_periods, max_periods=max_periods, data_asdf_directory=data_asdf_directory,
                              sync_asdf_directory=sync_asdf_directory, data_info_directory=data_info_directory)
        gcmtids_this_rank = get_gcmtids_this_rank(windows_directory)
        for each_gcmtid in gcmtids_this_rank:
            misfit_windows = kernel_used(each_gcmtid)
            save_pickle_event(misfit_windows, output_directory, each_gcmtid)

    main()  # pylint: disable=no-value-for-parameter
