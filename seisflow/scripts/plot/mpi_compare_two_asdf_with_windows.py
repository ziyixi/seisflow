"""
mpi_compare_two_asdf_with_windows.py: plot the whole trace waveform comparision pdfs.
"""
from glob import glob
from os.path import basename, join

import click
import numpy as np
from mpi4py import MPI

from ...setting import SURFACE_THRESHOLD
from ...utils.load_files import load_pickle
from .compare_two_asdf_with_windows import kernel

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
@click.option('--misfit_windows_dir', required=True, type=str)
@click.option('--obs_dir', required=True, type=str)
@click.option('--syn_dir', required=True, type=str)
@click.option('--output_dir', required=True, type=str)
@click.option('--info_dir', required=True, type=str)
@click.option('--azimuth_width', required=True, type=int)
@click.option('--waves_perpage', required=True, type=int)
@click.option('--snr', required=True, type=float)
@click.option('--cc', required=True, type=float)
@click.option('--deltat', required=True, type=float)
@click.option('--band', required=True, type=str)
@click.option('--plot_surface', is_flag=True, default=False)
def main(misfit_windows_dir, obs_dir, syn_dir, output_dir, info_dir, azimuth_width, waves_perpage, snr, cc, deltat, band, plot_surface):
    all_gcmtids = get_all_gcmtids_with_filtered(
        misfit_windows_dir, info_dir, plot_surface)
    gcmtids_this_rank = get_gcmtids_each_rank(all_gcmtids)
    min_time, max_time = map(int, band.split(","))
    for each_gcmtid in gcmtids_this_rank:
        # * build up parameters
        obs_asdf = join(
            obs_dir, f"{each_gcmtid}.preprocessed_{min_time}s_to_{max_time}s.h5")
        syn_asdf = join(
            syn_dir, f"{each_gcmtid}.preprocessed_{min_time}s_to_{max_time}s.h5")
        output_pdf = join(
            output_dir, f"{each_gcmtid}.preprocessed_{min_time}s_to_{max_time}s.pdf")
        misfit_windows_path = join(misfit_windows_dir, f"{each_gcmtid}.pkl")

        # * run the kernel
        kernel(obs_asdf, syn_asdf, azimuth_width, output_pdf, waves_perpage,
               info_dir, misfit_windows_path, snr, cc, deltat, False, plot_surface)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
