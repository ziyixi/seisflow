"""
line_search.py: do the line search for single event.
"""

from copy import copy

import click
import numpy as np

from ...setting import LINE_SEARCH_PERTURBATION_BOUNDARY
from ...tasks.line_search.line_search_structure import \
    calculate_weighted_misfit
from ...utils.asdf_io import VirAsdf
from ...utils.load_files import load_first_arrival_baz_evdp, load_pickle
from .update_model_stw105 import kernel as update_model_stw105_kernel


def get_perturbed_vir_sync(virasdf_raw, virasdf_perturbed, step_length):
    """
    get_perturbed_vir_sync: get perturbed virasdf given a step length.
    """
    virasdf_vir_perturbed = copy(virasdf_raw)
    for each_net_sta in virasdf_vir_perturbed.get_waveforms_list():
        st_raw = virasdf_raw.get_waveforms()[each_net_sta]["st"]
        st_per = virasdf_perturbed.get_waveforms()[each_net_sta]["st"]
        st_vir_perturbed = virasdf_vir_perturbed.get_waveforms()[
            each_net_sta]["st"]
        for index in range(len(st_vir_perturbed)):
            st_vir_perturbed[index].data = st_raw[index].data + \
                (step_length/LINE_SEARCH_PERTURBATION_BOUNDARY) * \
                (st_per[index].data-st_raw[index].data)
        virasdf_vir_perturbed.update_st(st_vir_perturbed, each_net_sta)
    return virasdf_vir_perturbed


def kernel(windows_path, data_asdf_path, sync_asdf_path_raw, sync_asdf_path_perturbed, data_info_directory, stations_path,
           search_range, search_step, lat_range, lon_range, kernel_path, model_path, output_path_d660, output_path_ppm):
    # * prepare to used parameters
    windows = load_pickle(windows_path)
    consider_surface = False
    data_virasdf_body = VirAsdf()
    data_virasdf_body.read_asdf(data_asdf_path)
    sync_virasdf_body_raw = VirAsdf()
    sync_virasdf_body_raw.read_asdf(sync_asdf_path_raw)
    sync_virasdf_body_perturbed = VirAsdf()
    sync_virasdf_body_perturbed.read_asdf(sync_asdf_path_perturbed)
    data_virasdf_surface = None
    sync_virasdf_surface = None
    first_arrival_zr, first_arrival_t, baz, _ = load_first_arrival_baz_evdp(
        data_info_directory)
    stations = np.loadtxt(stations_path, dtype=np.str)

    # * now we can do the line search
    search_min, search_max = map(float, search_range.split(","))
    search_list = np.arange(search_min, search_max+search_step, search_step)
    # * for each searching value, we calculate the misfit
    weighted_misfit_collection = []
    for each_search_step in search_list:
        sync_virasdf_body_per = get_perturbed_vir_sync(
            sync_virasdf_body_raw, sync_virasdf_body_perturbed, each_search_step)
        weighted_misfit_collection.append(calculate_weighted_misfit(windows, consider_surface, data_virasdf_body, sync_virasdf_body_per, data_virasdf_surface,
                                                                    sync_virasdf_surface, first_arrival_zr, first_arrival_t, baz, stations))

    # * print out the searching result and return the search length
    print("line search result:")
    for index, each_search_step in enumerate(search_list):
        print(
            f"step: {each_search_step}, misfit: {weighted_misfit_collection[index]}")
    min_weighted_misfit = np.min(weighted_misfit_collection)
    optimized_step_length = search_list[np.argmin(
        weighted_misfit_collection)]
    print("optimized value:")
    print(
        f"step: {optimized_step_length}, misfit: {min_weighted_misfit}")
    # * now we generate the optimized d660 file and the ppm file.
    update_model_stw105_kernel(lat_range, lon_range, kernel_path,
                               model_path, output_path_d660, output_path_ppm, optimized_step=optimized_step_length)


@click.command()
@click.option('--windows_path', required=True, type=str, help="the windows path")
@click.option('--data_asdf_path', required=True, type=str, help="the data asdf path")
@click.option('--sync_asdf_path_raw', required=True, type=str, help="the raw processed sync asdf path")
@click.option('--sync_asdf_path_perturbed', required=True, type=str, help="the perturbed processed sync asdf path")
@click.option('--data_info_directory', required=True, type=str, help="the data info directory")
@click.option('--stations_path', required=True, type=str, help="the stations path")
@click.option('--search_range', required=True, type=str, help="search_min,search_max")
@click.option('--search_step', required=True, type=float, help="the search step")
@click.option('--lat_range', required=True, type=str, help="latmin,latmax should be larger than the simulation region")
@click.option('--lon_range', required=True, type=str, help="lonmin,lonmax should be larger than the simulation region")
@click.option('--kernel_path', required=True, type=str, help="the generated kernel text file from julia")
@click.option('--model_path', required=True, type=str, help="the external mesh model file")
@click.option('--output_path_d660', required=True, type=str, help="the external mesh file output path")
@click.option('--output_path_ppm', required=True, type=str, help="the external ppm file output path")
def main(windows_path, data_asdf_path, sync_asdf_path_raw, sync_asdf_path_perturbed, data_info_directory, stations_path,
         search_range, search_step, lat_range, lon_range, kernel_path, model_path, output_path_d660, output_path_ppm):
    kernel(windows_path, data_asdf_path, sync_asdf_path_raw, sync_asdf_path_perturbed, data_info_directory, stations_path,
           search_range, search_step, lat_range, lon_range, kernel_path, model_path, output_path_d660, output_path_ppm)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
