"""
mpi_compare_phases_improvements.py: generate the phases comparision pdf files.
"""
from glob import glob
from os.path import basename, join

import click
import numpy as np
from mpi4py import MPI

from ...setting import SURFACE_THRESHOLD
from ...utils.load_files import load_pickle
from .compare_phases_improvements import kernel

comm = MPI.COMM_WORLD  # pylint: disable=c-extension-no-member
size = comm.Get_size()
rank = comm.Get_rank()


def get_all_gcmtids_with_filtered(misfit_windows_dir, info_dir, plot_surface):
    all_misfit_windows_path = sorted(glob(join(misfit_windows_dir, "*pkl")))
    all_gcmtids = [basename(item).split(".")[0]
                   for item in all_misfit_windows_path]
    evdp_dict = load_pickle(join(info_dir, "extra.evdp.pkl"))
    surface_result = {}
    for each_gcmtid in all_gcmtids:
        rep_net_sta = list(evdp_dict[each_gcmtid].keys())[0]
        depth = evdp_dict[each_gcmtid][rep_net_sta]
        if (depth <= SURFACE_THRESHOLD):
            surface_result[each_gcmtid] = True
        else:
            surface_result[each_gcmtid] = False
    if (plot_surface):
        final_gcmtids = []
        for each_gcmtid in all_gcmtids:
            if (surface_result[each_gcmtid]):
                final_gcmtids.append(each_gcmtid)
    else:
        final_gcmtids = all_gcmtids
    return final_gcmtids


def get_gcmtids_each_rank(all_gcmtids):
    return np.array_split(all_gcmtids, size)[rank]


@click.command()
@click.option('--obs_dir', required=True, type=str)
@click.option('--syn_dir1', required=True, type=str)
@click.option('--syn_dir2', required=True, type=str)
@click.option('--azimuth_width', required=True, type=int)
@click.option('--output_dir', required=True, type=str)
@click.option('--waves_perpage', required=True, type=int)
@click.option('--info_dir', required=True, type=str)
@click.option('--misfit_windows_dir1', required=True, type=str)
@click.option('--misfit_windows_dir2', required=True, type=str)
@click.option('--snr', required=True, type=float)
@click.option('--cc', required=True, type=float)
@click.option('--deltat', required=True, type=float)
@click.option('--band', required=True, type=str)
@click.option('--plot_surface', is_flag=True, default=False)
def main(obs_dir, syn_dir1, syn_dir2, azimuth_width, output_dir, waves_perpage, info_dir, misfit_windows_dir1,
         misfit_windows_dir2, snr, cc, deltat, band, plot_surface):
    all_gcmtids = get_all_gcmtids_with_filtered(
        misfit_windows_dir1, info_dir, plot_surface)
    gcmtids_this_rank = get_gcmtids_each_rank(all_gcmtids)
    min_time, max_time = map(int, band.split(","))
    for each_gcmtid in gcmtids_this_rank:
        # * build up parameters
        obs_asdf = join(
            obs_dir, f"{each_gcmtid}.preprocessed_{min_time}s_to_{max_time}s.h5")
        syn_asdf_1 = join(
            syn_dir1, f"{each_gcmtid}.preprocessed_{min_time}s_to_{max_time}s.h5")
        syn_asdf_2 = join(
            syn_dir2, f"{each_gcmtid}.preprocessed_{min_time}s_to_{max_time}s.h5")
        misfit_windows_path_1 = join(
            misfit_windows_dir1, f"{each_gcmtid}.pkl")
        misfit_windows_path_2 = join(
            misfit_windows_dir2, f"{each_gcmtid}.pkl")
        # * now we call the plotting script to plot
        kernel(obs_asdf, syn_asdf_1, syn_asdf_2, azimuth_width, output_dir, waves_perpage,
               info_dir, misfit_windows_path_1, misfit_windows_path_2, snr, cc, deltat, False, plot_surface)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
