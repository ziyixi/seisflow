"""
calculate_adjoint_source_zerolagcc_one_event: calculate the adjoint source for a single event.
"""
from collections import namedtuple

import numpy as np

from ...utils.load_files import get_stations_mapper
from .adjoint_source_each_window_zerolagcc import \
    calculate_adjoint_source_each_window
from .weight import (cal_category_weight, cal_cc_weight, cal_deltat_weight,
                     cal_geographical_weight, cal_snr_weight)

Weight = namedtuple(
    'Weight', ['snr', 'cc', 'deltat', 'geographical', 'category'])


def get_weights_for_all(misfit_windows, stations,  snr_threshold, cc_threshold, deltat_threshold, calculate_basic):
    """
    get_weights_for_all: calculate weights.
    """
    weights_for_all = {}
    # * firstly we update the weight of snr,cc,deltat
    for net_sta in misfit_windows:
        weights_for_all[net_sta] = {}
        for category in misfit_windows[net_sta]:
            weights_for_all[net_sta][category] = []
            for each_misfit_window in misfit_windows[net_sta][category].windows:
                wsnr = cal_snr_weight(each_misfit_window,
                                      snr_threshold[0], snr_threshold[1])
                wcc = cal_cc_weight(each_misfit_window,
                                    cc_threshold[0], cc_threshold[1])
                wdeltat = cal_deltat_weight(each_misfit_window,
                                            deltat_threshold[0], deltat_threshold[1])
                weights_for_all[net_sta][category].append(
                    Weight(wsnr, wcc, wdeltat, None, None))
    if(not calculate_basic):
        # * get the station list for the geographical weighting (remove all 0 cases)
        used_geographical_net_sta_list = []
        for net_sta in weights_for_all:
            status = False
            for category in weights_for_all[net_sta]:
                for each_weight in weights_for_all[net_sta][category]:
                    wsnr_cc_deltat = each_weight.snr * each_weight.cc * each_weight.deltat
                    if (wsnr_cc_deltat > 0):
                        status = True
            if (status):
                used_geographical_net_sta_list.append(net_sta)
        # build stations_mapper
        stations_mapper = get_stations_mapper(stations)
        # get geographical weighting and update
        geographical_weight_dict = cal_geographical_weight(
            stations_mapper, used_geographical_net_sta_list, list(weights_for_all.keys()))
        for net_sta in weights_for_all:
            for category in weights_for_all[net_sta]:
                for index, each_weight in enumerate(weights_for_all[net_sta][category]):
                    weights_for_all[net_sta][category][index] = each_weight._replace(
                        geographical=geographical_weight_dict[net_sta])

        # * get the number of items for each category
        # firstly we get all the category names
        rep_net_sta = list(weights_for_all.keys())[0]
        all_categories = list(weights_for_all[rep_net_sta].keys())
        # here we should weight based on number of windows but not the number of usable stations.
        number_each_category = {}
        for each_category in all_categories:
            number_each_category[each_category] = 0
            for net_sta in weights_for_all:
                for each_weight in weights_for_all[net_sta][each_category]:
                    # if this window is usable or not
                    wsnr_cc_deltat = each_weight.snr * each_weight.cc * each_weight.deltat
                    if (wsnr_cc_deltat > 0):
                        number_each_category[each_category] += 1
        # get category weighting and update
        # here we should weight based on number of windows but not the number of usable stations.
        weight_each_category = {}
        for each_category in number_each_category:
            weight_each_category[each_category] = cal_category_weight(
                number_each_category[each_category])
        for net_sta in weights_for_all:
            for category in weights_for_all[net_sta]:
                for index, each_weight in enumerate(weights_for_all[net_sta][category]):
                    weights_for_all[net_sta][category][index] = each_weight._replace(
                        category=weight_each_category[category])
    return weights_for_all


def calculate_adjoint_source_zerolagcc_one_event(misfit_windows, stations, raw_sync_virasdf, snr_threshold, cc_threshold, deltat_threshold, body_band, surface_band,
                                                 consider_surface, sync_virasdf_body, data_virasdf_body, sync_virasdf_surface, data_virasdf_surface):
    """
    calculate_adjoint_source_zerolagcc_one_event
        + misfit_windows: misfit_windows[net_sta][category_name] as Windows_collection
        + stations: stations in Specfem3D format
        + raw_virsync_asdf: raw vir sync asdf 
        + snr_threshold: (snr_value1,snr_value2)
        + cc_threshold: (cc_value1,cc_value2)
        + deltat_threshold: (deltat_value1,deltat_value2)
        + body_band: body filter band in seconds (mintime,maxtime)
        + surface_band: surface filter band in seconds (mintime,maxtime)
        + consider_surface: if consider the surface waves
        + sync_virasdf_body: sync vir asdf for the body wave
        + data_virasdf_body: data vir asdf for the body wave
        + sync_virasdf_surface: sync vir asdf for the surface wave
        + data_virasdf_surface: data vir asdf for the surface wave
    output:
        a dict result, with result[net_sta] as a numpy array of shape (3, len(trace.data))
    """
    # * get weight for all the windows with the same structure as misfit_windows
    weights_for_all = get_weights_for_all(
        misfit_windows, stations, snr_threshold, cc_threshold, deltat_threshold, False)

    # * get adjoint sources (weighted)
    # weights_for_all has the same structure with misfit_windows
    adjoint_source_zerolagcc = {}
    weight_normalize_factor = 0
    # get adjoint_source_length
    rep_net_sta = list(weights_for_all.keys())[0]
    rep_sync_wg = raw_sync_virasdf.get_waveforms()[rep_net_sta]
    rep_trace = rep_sync_wg["st"][0]
    adjoint_source_length = len(rep_trace.data)

    # ! fix a bug, when the event time of raw_sync_virasdf and sync_virasdf_body are different, we have to consider raw_sync_asdf_trace_adjust_time
    true_event_time = sync_virasdf_body.get_events()[0].preferred_origin().time
    raw_event_time = raw_sync_virasdf.get_events()[0].preferred_origin().time
    raw_sync_asdf_trace_adjust_time = true_event_time-raw_event_time

    for net_sta in weights_for_all:
        # in the order of E,T,Z
        adjoint_source_zerolagcc[net_sta] = np.zeros(
            (3, adjoint_source_length))
        for category in weights_for_all[net_sta]:
            if (category in ["z", "r", "t"]):
                data_asdf = data_virasdf_body
                sync_asdf = sync_virasdf_body
                mintime, maxtime = body_band
            elif (category in ["surface_z", "surface_r", "surface_t"]):
                data_asdf = data_virasdf_surface
                sync_asdf = sync_virasdf_surface
                mintime, maxtime = surface_band
            data_wg = data_asdf.get_waveforms()[net_sta]
            sync_wg = sync_asdf.get_waveforms()[net_sta]
            raw_sync_wg = raw_sync_virasdf.get_waveforms()[net_sta]
            for index, each_misfit_window in enumerate(misfit_windows[net_sta][category].windows):
                each_weight = weights_for_all[net_sta][category][index]
                wcc = each_weight.cc
                wdeltat = each_weight.deltat
                wsnr = each_weight.snr
                wcategory = each_weight.category
                # wgeographical = each_weight.geographical
                # ! sync test, set geo weighting as 1
                wgeographical = 1.0
                weight_normalize_factor += wcc * wdeltat * wsnr * wcategory * wgeographical

                # calculate each_adjoint_trace
                component = each_misfit_window.component
                sync_asdf_trace = sync_wg["st"].select(
                    component=component)[0]
                raw_sync_asdf_trace = raw_sync_wg["st"].select(component=component)[
                    0]
                data_asdf_trace = data_wg["st"].select(
                    component=component)[0]
                each_adjoint_trace = calculate_adjoint_source_each_window(
                    each_misfit_window, raw_sync_asdf_trace, sync_asdf_trace, data_asdf_trace, mintime, maxtime, raw_sync_asdf_trace_adjust_time=raw_sync_asdf_trace_adjust_time)

                if (component == "Z"):
                    adjoint_source_zerolagcc[net_sta][2, :] += each_adjoint_trace.data * \
                        wcc * wdeltat * wsnr * wcategory * wgeographical
                elif (component == "R"):
                    r_adjoint_source = each_adjoint_trace.data * \
                        wcc * wdeltat * wsnr * wcategory * wgeographical
                    theta = (each_misfit_window.baz - 180) % 360
                    adjoint_source_zerolagcc[net_sta][0, :] += r_adjoint_source * \
                        np.sin(np.deg2rad(theta))
                    adjoint_source_zerolagcc[net_sta][1, :] += r_adjoint_source * \
                        np.cos(np.deg2rad(theta))
                elif (component == "T"):
                    t_adjoint_source = each_adjoint_trace.data*wcc * \
                        wdeltat * wsnr * wcategory * wgeographical
                    theta = (each_misfit_window.baz - 90) % 360
                    adjoint_source_zerolagcc[net_sta][0, :] += t_adjoint_source * \
                        np.sin(np.deg2rad(theta))
                    adjoint_source_zerolagcc[net_sta][1, :] += t_adjoint_source * \
                        np.cos(np.deg2rad(theta))
    # normalize the adjoint source
    for net_sta in adjoint_source_zerolagcc:
        adjoint_source_zerolagcc[net_sta] /= weight_normalize_factor
        # check if the adjoint source for this station is nan
        # now we have only seen the station HB.ZHX has some problem
        status = np.isfinite(adjoint_source_zerolagcc[net_sta]).all()
        if (not status):
            # there is inf/-inf/nan in the adjoint source
            adjoint_source_zerolagcc[net_sta][:] = 0.0

    return adjoint_source_zerolagcc
