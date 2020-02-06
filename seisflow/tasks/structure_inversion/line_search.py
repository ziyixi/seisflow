"""
line_search.py: perform the line search for the structure inversion.
"""
from copy import copy

import numpy as np

from ...utils.setting import LINE_SEARCH_PERTURBATION, SNR_THRESHOLD, CC_THRESHOLD, DELTAT_THRESHOLD
from ..adjoint.calculate_adjoint_source_zerolag_cc_multiple_events import \
    get_weights_for_all
from ..windows.calculate_misfit_windows import calculate_misfit_windows
from mpi4py import MPI


class Store():
    weight = None


comm = MPI.COMM_WORLD  # pylint: disable=c-extension-no-member
size = comm.Get_size()
rank = comm.Get_rank()


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
                (step_length/LINE_SEARCH_PERTURBATION) * \
                (st_per[index].data-st_raw[index].data)
        virasdf_vir_perturbed.update_st(st_vir_perturbed, each_net_sta)
    return virasdf_vir_perturbed


def calculate_weighted_misfit(windows, consider_surface, data_virasdf_body, sync_virasdf_body, data_virasdf_surface, sync_virasdf_surface, first_arrival_zr, first_arrival_t, baz,
                              stations):
    """
    calculate_weighted_misfit: calculate the weighted misfit for a singlue event.
    """
    # * firstly we calculate the misfit windows
    misfit_windows = calculate_misfit_windows(windows, consider_surface, data_virasdf_body, sync_virasdf_body,
                                              data_virasdf_surface, sync_virasdf_surface, first_arrival_zr, first_arrival_t, baz)
    # * as we have to sum all the weighted misfits and normalization factors for all the event, we have to return two values
    if (Store.weight == None):
        if(rank == 1):
            print("start to calculate misfit.")
        snr_threshold = tuple(map(float, SNR_THRESHOLD.split(",")))
        cc_threshold = tuple(map(float, CC_THRESHOLD.split(",")))
        deltat_threshold = tuple(map(float, DELTAT_THRESHOLD.split(",")))
        # as we want to calculate the geographical and categorical weighting
        calculate_basic = False
        # note get_weights_for_all is from weighting calculation for multiple events, thus use mpi here. [parallel]
        Store.weight = get_weights_for_all(misfit_windows, stations, snr_threshold,
                                           cc_threshold, deltat_threshold, calculate_basic)
    summation_weighted_misfit = 0
    summation_weighting_factor = 0
    for net_sta in misfit_windows:
        for category in misfit_windows[net_sta]:
            for index in range(len(misfit_windows[net_sta][category].windows)):
                each_misfit_window = misfit_windows[net_sta][category].windows[index]
                # unsubscriptable-object
                each_weight = Store.weight[net_sta][category][index]  # pylint: disable=unsubscriptable-object
                # get weighted similarity for this window
                wcc = each_weight.cc
                wdeltat = each_weight.deltat
                wsnr = each_weight.snr
                wcategory = each_weight.category
                wgeographical = each_weight.geographical
                wtotal = wcc * wdeltat * wsnr * wcategory * wgeographical
                each_similarity = each_misfit_window.similarity
                if (np.isnan(each_similarity)):
                    # don't count this window
                    wtotal = 0
                    each_similarity = 0
                summation_weighted_misfit += (1-each_similarity)*wtotal
                summation_weighting_factor += wtotal
    # * here we should gather the weighted misfit and weighting factor for all the events
    summation_weighted_misfit_all_events = mpi_gather_summation_weighted_misfit(
        summation_weighted_misfit)
    summation_weighting_factor_all_events = mpi_gather_summation_weighting_factor(
        summation_weighting_factor)
    final_weighted_misfit = summation_weighted_misfit_all_events / \
        summation_weighting_factor_all_events
    return final_weighted_misfit


def mpi_gather_summation_weighted_misfit(summation_weighted_misfit):
    """
    mpi_gather_summation_weighted_misfit: gather weighted misfit for all the events
    """
    summation_weighted_misfit_all_events_list = comm.allgather(
        summation_weighted_misfit)
    comm.Barrier()
    summation_weighted_misfit_all_events = np.sum(
        summation_weighted_misfit_all_events_list)
    return summation_weighted_misfit_all_events


def mpi_gather_summation_weighting_factor(summation_weighting_factor):
    """
    mpi_gather_summation_weighting_factor: gather weighting factor for all the events:
    """
    summation_weighting_factor_all_events_list = comm.allgather(
        summation_weighting_factor)
    comm.Barrier()
    summation_weighting_factor_all_events = np.sum(
        summation_weighting_factor_all_events_list)
    return summation_weighting_factor_all_events


def line_search(data_virasdf_body, data_virasdf_surface, sync_virasdf_body_raw, sync_virasdf_surface_raw, sync_virasdf_body_perturbed, sync_virasdf_surface_perturbed,
                windows, consider_surface, first_arrival_zr, first_arrival_t, baz, stations,
                search_range, search_step_space):
    """
    line_search: perform the line search to find the optimized step length for model updating.
    """
    # ! note here the meaning of raw is processed raw, not the same as the source inversion.
    # the processing part should be done before doing the line search (as we can directly process the non-green waveforms)
    search_step_lengths = np.arange(
        search_range[0], search_range[1] + search_step_space, search_step_space)
    weighted_misfit_collection = []
    for each_search_step_length in search_step_lengths:
        # * for the given step length, firstly we have to calculate the perturbed sync
        sync_virasdf_body = get_perturbed_vir_sync(
            sync_virasdf_body_raw, sync_virasdf_body_perturbed, each_search_step_length)
        sync_virasdf_surface = get_perturbed_vir_sync(
            sync_virasdf_surface_raw, sync_virasdf_surface_perturbed, each_search_step_length)
        # * for each search step length, we calculate the summed weighted misfit value
        weighted_misfit_collection.append(calculate_weighted_misfit(windows, consider_surface, data_virasdf_body, sync_virasdf_body, data_virasdf_surface,
                                                                    sync_virasdf_surface, first_arrival_zr, first_arrival_t, baz, stations))
        # * we want to find the min value of the weighted misfit and corresponding step length
    min_weighted_misfit = np.min(weighted_misfit_collection)
    optimized_step_length = search_step_lengths[np.argmin(
        weighted_misfit_collection)]
    return min_weighted_misfit, optimized_step_length
