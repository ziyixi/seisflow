"""
source_line_search.py: give the waveforms from green func and the data, get optimized parameters and perturbed new CMTSOLUTION.
"""
from glob import glob
from os.path import basename, join

import numpy as np
import obspy
from mpi4py import MPI

from ..tasks.source.line_search import line_search
from ..tasks.source.make_perturbed_cmtsolution import add_src_frechet
from ..utils.asdf_io import VirAsdf
from ..utils.get_path import get_data_asdf_fnames
from ..utils.load_files import load_first_arrival_baz_evdp, load_pickle_event

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


def load_src_frechet(src_frechet_path):
    return np.loadtxt(src_frechet_path)


def load_cmtsolution(cmtsolution_path):
    cmtsolution = obspy.read_events(cmtsolution_path)[0]
    return cmtsolution


# * ================================================================
# we have the following group of parameters
# * paths
#   green_raw_asdf_directory
#   green_perturbed_asdf_directory
#   data_asdf_directory
#   windows directory
#   data_info_directory
#   stations_path
#   src_frechet_directory
#   cmtsolution_directory
#   output_newcmtsolution_directory
# * parameters used to process the sync, also used to get data paths
#   body_band
#   surface_band
#   taper_tmin_tmax
# * parameters used in the optimization
#   random_state
#   alpha_range
#   t0_range
#   tau_range (should only be changed in the last step or not to change)
#   tau_raw
#   t0_raw
#   init_points
#   n_iter
# * parameters used to set weight and other parameters
#   snr_threshold
#   cc_threshold
#   deltat_threshold
#   max_dxs_ratio
# * ================================================================
def get_paths(green_raw_asdf_directory, green_perturbed_asdf_directory, data_asdf_directory, windows_directory, data_info_directory, stations_path,
              src_frechet_directory, cmtsolution_directory, output_newcmtsolution_directory, body_band, surface_band):
    """
    get_paths: handle all the paths related part
    """
    # * get all the gcmtids used, use files in cmtsolution_directory as the reference
    all_files_cmtsolution_directory = glob(join(cmtsolution_directory, "*"))
    all_gcmtids = [basename(item) for item in all_files_cmtsolution_directory]
    all_gcmtids_this_rank = np.array_split(all_gcmtids, size)[rank]
    # * get green_raw_asdf_directory for this rank, convert to virasdf
    all_files_green_raw_asdf_directory_this_rank = [
        join(green_raw_asdf_directory, gcmtid+".h5") for gcmtid in all_gcmtids_this_rank]
    all_virasdf_green_raw_this_rank = []
    for each_file in all_files_green_raw_asdf_directory_this_rank:
        virasdf = VirAsdf()
        virasdf.read_asdf(each_file)
        all_virasdf_green_raw_this_rank.append(virasdf)
    # * get green_perturbed_asdf_directory for this rank, convert to virasdf
    all_files_green_perturbed_asdf_directory_this_rank = [
        join(green_perturbed_asdf_directory, gcmtid+".h5") for gcmtid in all_gcmtids_this_rank]
    all_virasdf_green_perturbed_this_rank = []
    for each_file in all_files_green_perturbed_asdf_directory_this_rank:
        virasdf = VirAsdf()
        virasdf.read_asdf(each_file)
        all_virasdf_green_perturbed_this_rank.append(virasdf)
    # * data_asdf_body_path, data_asdf_surface_path and convert them to virasdf
    all_virasdf_data_body_this_rank = []
    all_virasdf_data_surface_this_rank = []
    min_periods = f"{body_band[0]},{surface_band[0]}"
    max_periods = f"{body_band[1]},{surface_band[1]}"
    for each_gcmtid in all_gcmtids_this_rank:
        data_asdf_body_path, data_asdf_surface_path = get_data_asdf_fnames(
            each_gcmtid, min_periods, max_periods, data_asdf_directory)
        data_virasdf_body = VirAsdf()
        data_virasdf_body.read_asdf(data_asdf_body_path)
        data_virasdf_surface = VirAsdf()
        data_virasdf_surface.read_asdf(data_asdf_surface_path)
        all_virasdf_data_body_this_rank.append(data_virasdf_body)
        all_virasdf_data_surface_this_rank.append(data_virasdf_surface)
    # * load windows
    all_windows_this_rank = []
    for each_gcmtid in all_gcmtids_this_rank:
        all_windows_this_rank.append(
            load_pickle_event(windows_directory, each_gcmtid))
    # * load first arrival, baz and evdp
    first_arrival_zr, first_arrival_t, baz, evdp = load_first_arrival_baz_evdp(
        data_info_directory)
    # * load stations
    stations = np.loadtxt(stations_path)
    # * get src_frechet used in this rank
    all_files_src_frechet_directory_this_rank = [
        join(src_frechet_directory, item) for item in all_gcmtids_this_rank]
    # * get cmtsolution used in this rank
    all_files_cmtsolution_directory_this_rank = [
        join(cmtsolution_directory, item) for item in all_gcmtids_this_rank]
    # * get output_newcmtsolution_directory path used in this rank
    all_files_output_newcmtsolution_directory_this_rank = [
        join(output_newcmtsolution_directory, item) for item in all_gcmtids_this_rank]
    # * return the results
    return all_virasdf_green_raw_this_rank, all_virasdf_green_perturbed_this_rank,
    all_virasdf_data_body_this_rank, all_virasdf_data_surface_this_rank, all_windows_this_rank,
    first_arrival_zr, first_arrival_t, baz, evdp, stations, all_files_src_frechet_directory_this_rank,
    all_files_cmtsolution_directory_this_rank, all_files_output_newcmtsolution_directory_this_rank


def main(green_raw_asdf_directory, green_perturbed_asdf_directory, data_asdf_directory, windows_directory, data_info_directory, stations_path,
         src_frechet_directory, cmtsolution_directory, output_newcmtsolution_directory,
         body_band, surface_band, taper_tmin_tmax,
         random_state, alpha_range, t0_range, tau_range, init_points, n_iter,
         snr_threshold, cc_threshold, deltat_threshold, max_dxs_ratio):
    # handle paths
    all_virasdf_green_raw_this_rank, all_virasdf_green_perturbed_this_rank, all_virasdf_data_body_this_rank, \
        all_virasdf_data_surface_this_rank, all_windows_this_rank, first_arrival_zr, first_arrival_t, baz, evdp, \
        stations, all_files_src_frechet_directory_this_rank, all_files_cmtsolution_directory_this_rank, \
        all_files_output_newcmtsolution_directory_this_rank = get_paths(green_raw_asdf_directory, green_perturbed_asdf_directory,
                                                                        data_asdf_directory, windows_directory, data_info_directory, stations_path,
                                                                        src_frechet_directory, cmtsolution_directory, output_newcmtsolution_directory, body_band, surface_band)
    # convert input parameters.
