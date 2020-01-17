"""
mpi_source_line_search.py: give the waveforms from green func and the data, get optimized parameters and perturbed new CMTSOLUTION.
"""
from glob import glob
from os.path import basename, join

import click
import numpy as np
import obspy
from mpi4py import MPI

from ..tasks.source.line_search import line_search
from ..tasks.source.make_perturbed_cmtsolution import add_src_frechet
from ..utils.asdf_io import VirAsdf
from ..utils.get_path import get_data_asdf_fnames
from ..utils.load_files import load_first_arrival_baz_evdp, load_pickle_event
from ..utils.save_files import save_cmtsolution
from ..utils.setting import (CC_THRESHOLD, DELTAT_THRESHOLD, INIT_POINTS,
                             MAX_DXS_RATIO, N_ITER, RANDOM_STATE,
                             SNR_THRESHOLD, SURFACE_THRESHOLD)

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
# * test is passed for this script on 01/07/2020


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
    all_files_cmtsolution_directory = sorted(
        glob(join(cmtsolution_directory, "*")))
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
    body_band = tuple(map(float, body_band.split(",")))
    surface_band = tuple(map(float, surface_band.split(",")))
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
    stations = np.loadtxt(stations_path, dtype=np.str)
    # * get src_frechet used in this rank
    all_files_src_frechet_directory_this_rank = [
        join(src_frechet_directory, item) for item in all_gcmtids_this_rank]
    # * get cmtsolution used in this rank
    all_files_cmtsolution_directory_this_rank = [
        join(cmtsolution_directory, item) for item in all_gcmtids_this_rank]
    # * get output_newcmtsolution_directory path used in this rank
    all_files_output_newcmtsolution_directory_this_rank = [
        join(output_newcmtsolution_directory, item) for item in all_gcmtids_this_rank]
    # * return the results, note we have to return tuple or it will only return the values in the first line (as C++)
    return (all_virasdf_green_raw_this_rank, all_virasdf_green_perturbed_this_rank,
            all_virasdf_data_body_this_rank, all_virasdf_data_surface_this_rank, all_windows_this_rank,
            first_arrival_zr, first_arrival_t, baz, evdp, stations, all_files_src_frechet_directory_this_rank,
            all_files_cmtsolution_directory_this_rank, all_files_output_newcmtsolution_directory_this_rank
            )


def convert_input_parameters(alpha_range, t0_range,
                             tau_range, body_band, surface_band, snr_threshold,
                             cc_threshold, deltat_threshold):
    """
    convert_input_parameters: convert the input form of some parameters to be used in the script
    """
    alpha_range = tuple(map(float, alpha_range.split(",")))
    t0_range = tuple(map(float, t0_range.split(",")))
    if (tau_range == "fixed"):
        pass
    else:
        tau_range = tuple(map(float, tau_range.split(",")))
    body_band = tuple(map(float, body_band.split(",")))
    surface_band = tuple(map(float, surface_band.split(",")))
    snr_threshold = tuple(map(float, snr_threshold.split(",")))
    cc_threshold = tuple(map(float, cc_threshold.split(",")))
    deltat_threshold = tuple(map(float, deltat_threshold.split(",")))
    return alpha_range, t0_range, tau_range, body_band, surface_band, snr_threshold, cc_threshold, deltat_threshold


@click.command()
@click.option('--green_raw_asdf_directory', required=True, type=str, help="the directory of the asdf files from the raw green function cmtsolution")
@click.option('--green_perturbed_asdf_directory', required=True, type=str, help="the directory of the asdf files from the perturbed green function cmtsolution")
@click.option('--data_asdf_directory', required=True, type=str, help="the directory of the processed data asdf files")
@click.option('--windows_directory', required=True, type=str, help="the directory of the windows files in the pickle format")
@click.option('--data_info_directory', required=True, type=str, help="the data info directory")
@click.option('--stations_path', required=True, type=str, help="the stations file path in the specfem3D format")
@click.option('--src_frechet_directory', required=True, type=str, help="the src frechet files directory")
@click.option('--cmtsolution_directory', required=True, type=str, help="the cmtsolution directory")
@click.option('--output_newcmtsolution_directory', required=True, type=str, help="the new cmtsolution directory to output for the next iterations")
@click.option('--body_band', required=True, type=str, help="body wave filtering range: min_period,max_period")
@click.option('--surface_band', required=True, type=str, help="surface wave filtering range: min_period,max_period")
@click.option('--taper_tmin_tmax', required=True, type=str, help="frequency taper range: min_period,max_period")
@click.option('--alpha_range', required=True, type=str, help="alpha range: min_alpha,max_alpha")
@click.option('--t0_range', required=True, type=str, help="t0 range: min_t0,max_t0")
@click.option('--tau_range', required=True, type=str, help="tau range: min_tau,max_tau")
def main(green_raw_asdf_directory, green_perturbed_asdf_directory, data_asdf_directory, windows_directory, data_info_directory, stations_path,
         src_frechet_directory, cmtsolution_directory, output_newcmtsolution_directory,
         body_band, surface_band, taper_tmin_tmax,
         alpha_range, t0_range, tau_range):
    # handle paths
    all_virasdf_green_raw_this_rank, all_virasdf_green_perturbed_this_rank, all_virasdf_data_body_this_rank, \
        all_virasdf_data_surface_this_rank, all_windows_this_rank, first_arrival_zr, first_arrival_t, baz, evdp, \
        stations, all_files_src_frechet_directory_this_rank, all_files_cmtsolution_directory_this_rank, \
        all_files_output_newcmtsolution_directory_this_rank = get_paths(green_raw_asdf_directory, green_perturbed_asdf_directory,
                                                                        data_asdf_directory, windows_directory, data_info_directory, stations_path,
                                                                        src_frechet_directory, cmtsolution_directory, output_newcmtsolution_directory, body_band, surface_band)
    # convert input parameters.
    snr_threshold = SNR_THRESHOLD
    cc_threshold = CC_THRESHOLD
    deltat_threshold = DELTAT_THRESHOLD
    alpha_range, t0_range, tau_range, body_band, surface_band, snr_threshold, \
        cc_threshold, deltat_threshold = convert_input_parameters(alpha_range, t0_range,
                                                                  tau_range, body_band, surface_band, snr_threshold,
                                                                  cc_threshold, deltat_threshold)
    # call the line search script
    each_virasdf_combined = zip(
        all_virasdf_green_raw_this_rank, all_virasdf_green_perturbed_this_rank, all_virasdf_data_body_this_rank, all_virasdf_data_surface_this_rank, all_windows_this_rank,
        all_files_src_frechet_directory_this_rank, all_files_cmtsolution_directory_this_rank, all_files_output_newcmtsolution_directory_this_rank)
    for (each_virasdf_green_raw, each_virasdf_green_perturbed, each_virasdf_data_body,
         each_virasdf_data_surface, each_windows, each_src_frechet_path, each_cmtsolution_path, each_output_path) in each_virasdf_combined:
        # call the function of the line search
        # notice here the unit of depth is m, so we should divide it by 1000
        depth = each_virasdf_data_body.get_events(
        )[0].preferred_origin().depth / 1000.0
        gcmtid = each_windows.gcmtid
        if (depth <= SURFACE_THRESHOLD):
            consider_surface = True
        else:
            consider_surface = False

        each_cmtsolution = load_cmtsolution(each_cmtsolution_path)
        t0_raw = 0.0
        tau_raw = each_cmtsolution.focal_mechanisms[0].moment_tensor.source_time_function.duration / 2
        init_points = INIT_POINTS
        n_iter = N_ITER
        weighted_similarity_raw, optimizer_max = line_search(alpha_range, t0_range, tau_range, each_virasdf_green_raw, each_virasdf_green_perturbed, each_virasdf_data_body, each_virasdf_data_surface,
                                                             each_windows, consider_surface, body_band, surface_band, taper_tmin_tmax, first_arrival_zr, first_arrival_t, baz,
                                                             RANDOM_STATE, tau_raw, t0_raw, stations,  snr_threshold, cc_threshold, deltat_threshold, init_points=init_points, n_iter=n_iter)
        # print informations
        print(f"gcmtid: {gcmtid}")
        print(f"raw accuracy: {weighted_similarity_raw:.3f}")
        print(f"new accuracy: {optimizer_max['target']:.3f}")
        print(f"parameters combination:")
        print(
            f"alpha: {optimizer_max['params']['alpha']:.3f}; tau: {optimizer_max['params']['tau']:.3f}; t0: {optimizer_max['params']['t0']:.3f}")
        # generate the new cmtsolution
        each_src_frechet = np.loadtxt(each_src_frechet_path)
        max_dxs_ratio = MAX_DXS_RATIO
        optimized_dxs_ratio = -1*max_dxs_ratio*optimizer_max['params']['alpha']
        new_cmtsolution = add_src_frechet(
            each_src_frechet, each_cmtsolution, optimized_dxs_ratio)
        # update tau and event_time
        # ! a possible bug here when optimizer_max['params']['tau']==0: the half duration will be 1 in obspy.
        new_cmtsolution.focal_mechanisms[0].moment_tensor.source_time_function.duration = optimizer_max['params']['tau'] * 2
        new_cmtsolution.preferred_origin(
        ).time += optimizer_max['params']['t0']
        # write
        save_cmtsolution(each_output_path, new_cmtsolution)


if __name__ == "__main__":
    main()
