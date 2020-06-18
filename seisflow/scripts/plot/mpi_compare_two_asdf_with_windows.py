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


def main()
