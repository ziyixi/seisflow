"""
slurm_calculate_misfit_windows.py: calculate misfit windows
"""
from ..tasks.windows.calculate_misfit_windows import calculate_misfit_windows
from glob import glob
from os.path import join, basename
import pickle
from functools import partial
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


def load_pickle(pickle_path):
    with open(pickle_path, "rb") as f:
        data = pickle.load(f)
    return data


def load_windows(gcmtid, windows_directory):
    fname = join(windows_directory, f"{gcmtid}.pkl")
    windows = load_pickle(fname)
    return windows


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


def load_first_arrival_baz_evdp(data_info_directory):
    first_arrival_zr_path = join(
        data_info_directory, "traveltime.P.pkl")
    first_arrival_t_path = join(
        data_info_directory, "traveltime.S.pkl")
    baz_path = join(
        data_info_directory, "extra.baz.pkl")
    evdp_path = join(
        data_info_directory, "extra.evdp.pkl")
    first_arrival_zr = load_pickle(first_arrival_zr_path)
    first_arrival_t = load_pickle(first_arrival_t_path)
    baz = load_pickle(baz_path)
    evdp = load_pickle(evdp_path)
    return first_arrival_zr, first_arrival_t, baz, evdp


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


def save_misfit_windows(misfit_windows, output_dir, used_gcmtid):
    output_path = join(output_dir, f"{used_gcmtid}.pkl")
    with open(output_path, "wb") as f:
        pickle.dump(misfit_windows, f)


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
            save_misfit_windows(misfit_windows, output_directory, each_gcmtid)

    main()
