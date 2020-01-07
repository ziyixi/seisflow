"""
line_search.py: using BayesianOptimization to find the step length and the source parameters at the same time.
"""
from collections import namedtuple
from functools import partial

import numpy as np
from bayes_opt import BayesianOptimization

from ...utils.asdf_io import VirAsdf
from ..adjoint.calculate_adjoint_source_one_event import get_weights_for_all
from ..process.convolve_src_func import conv_sf_and_st, source_time_func
from ..process.process_sync_single_st import (ahead_process_green_function_st,
                                              post_process_green_function_st)
from ..windows.calculate_misfit_windows import calculate_misfit_windows

Weight = namedtuple(
    'Weight', ['snr', 'cc', 'deltat', 'geographical', 'category'])


def prepare_green_func(green_virasdf, body_band, surface_band, taper_tmin_tmax, consider_surface):
    """
    prepare_green_func: do the preprocesing for all the streams in virasdf
    """
    body_virasdf = green_virasdf.copy()
    if(consider_surface):
        surface_virasdf = green_virasdf.copy()
    else:
        surface_virasdf = None
    for each_net_sta in green_virasdf.get_waveforms_list():
        st_body = body_virasdf[each_net_sta]["st"]
        st_body = ahead_process_green_function_st(
            st_body, taper_tmin_tmax, body_band[0], body_band[1])
        body_virasdf.update_st(st_body, each_net_sta)
        if (consider_surface):
            st_surface = surface_virasdf[each_net_sta]["st"]
            st_surface = ahead_process_green_function_st(
                st_surface, taper_tmin_tmax, surface_band[0], surface_band[1])
            surface_virasdf.update_st(st_surface, each_net_sta)
    return body_virasdf, surface_virasdf


def get_perturbed_gf(green_virasdf1, green_virasdf2, green_virasdf_raw, alpha):
    """
    get_perturbed_gf: get asdf as waveforms=asdf1-asdf2 (usually zero-pos)
    """
    green_virasdf_perturbed = green_virasdf_raw.copy()
    for each_net_sta in green_virasdf_perturbed:
        st1 = green_virasdf1.get_waveforms()[each_net_sta]["st"]
        st2 = green_virasdf2.get_waveforms()[each_net_sta]["st"]
        st_perturbed = green_virasdf_perturbed.get_waveforms()[
            each_net_sta]["st"]
        for index in range(len(st_perturbed)):
            st_perturbed[index].data = st_perturbed[index].data + \
                alpha * (st1[index].data - st2[index].data)
        green_virasdf_perturbed.update_st(st_perturbed, each_net_sta)
    return green_virasdf_perturbed


def conv_and_postprocess(green_virasdf, sf, t0, waveform_length, sampling_rate):
    """
    conv_and_postprocess: do the convolution and postprocessing for all the sts in virasdf.
    """
    # we should not change green_virasdf, so make a copy of it.
    conv_virasdf = green_virasdf.copy()
    # firstly, we do the convolution
    for each_net_sta in conv_virasdf.get_waveforms_list():
        st = conv_virasdf.get_waveforms()[each_net_sta]["st"]
        st = conv_sf_and_st(st, sf, t0)
        conv_virasdf.update_st(st, each_net_sta)
    # and we can do the postprocessing
    event = green_virasdf.get_events()[0]
    origin = event.preferred_origin() or event.origins[0]
    event_time = origin.time
    for each_net_sta in conv_virasdf.get_waveforms_list():
        st = conv_virasdf.get_waveforms()[each_net_sta]["st"]
        st = post_process_green_function_st(
            st, event_time, waveform_length, sampling_rate)
        conv_virasdf.update_st(st, each_net_sta)
    return conv_virasdf


def forward_misfit_windows(alpha, t0, tau, body_green_virasdf1, body_green_virasdf2, body_green_virasdf_raw, surface_green_virasdf1,
                           surface_green_virasdf2, surface_green_virasdf_raw, data_virasdf_body, data_virasdf_surface,
                           windows, waveform_length, sampling_rate, consider_surface, first_arrival_zr, first_arrival_t, baz):
    """
    forward_misfit_windows: get misfit_windows in the forward simulation
    """
    # get perturbed virasdf
    body_green_virasdf_perturbed = get_perturbed_gf(
        body_green_virasdf1, body_green_virasdf2, body_green_virasdf_raw, alpha)
    if (consider_surface):
        surface_green_virasdf_perturbed = get_perturbed_gf(
            surface_green_virasdf1, surface_green_virasdf2, surface_green_virasdf_raw, alpha)
    else:
        surface_green_virasdf_perturbed = None
    # get sf
    rep_net_sta = body_green_virasdf_perturbed.get_waveforms_list()[0]
    dt = body_green_virasdf_perturbed.get_waveforms()[
        rep_net_sta]["st"].stats.delta
    sf = source_time_func(tau, dt)
    # conv and post process
    body_conv_virasdf_perturbed = conv_and_postprocess(
        body_green_virasdf_perturbed, sf, t0, waveform_length, sampling_rate)
    if (consider_surface):
        surface_conv_virasdf_perturbed = conv_and_postprocess(
            surface_green_virasdf_perturbed, sf, t0, waveform_length, sampling_rate)
    # calculate to get misfit windows
    misfit_windows = calculate_misfit_windows(windows, consider_surface, data_virasdf_body, body_conv_virasdf_perturbed,
                                              data_virasdf_surface, surface_conv_virasdf_perturbed, first_arrival_zr, first_arrival_t, baz)
    return misfit_windows


def forward_kernel(alpha, t0, tau, body_green_virasdf1, body_green_virasdf2, body_green_virasdf_raw, surface_green_virasdf1,
                   surface_green_virasdf2, surface_green_virasdf_raw, data_virasdf_body, data_virasdf_surface,
                   weights, windows, waveform_length, sampling_rate, consider_surface, first_arrival_zr, first_arrival_t, baz):
    """
    forward_kernel: the kernel function for the forward function used in the optimization
    """
    misfit_windows = forward_misfit_windows(alpha, t0, tau, body_green_virasdf1, body_green_virasdf2, body_green_virasdf_raw, surface_green_virasdf1,
                                            surface_green_virasdf2, surface_green_virasdf_raw, data_virasdf_body, data_virasdf_surface,
                                            windows, waveform_length, sampling_rate, consider_surface, first_arrival_zr, first_arrival_t, baz)
    weighted_similarity = get_weighted_similarity(misfit_windows, weights)
    return weighted_similarity


def get_weighted_similarity(misfit_windows, weights):
    """
    get_weighted_similarity: calculate the weighted misfit, return the weighted similarity
    """
    weight_normalize_factor = 0
    total_misfit = 0
    for net_sta in misfit_windows:
        for category in misfit_windows[net_sta]:
            for index in range(len(misfit_windows[net_sta][category].windows)):
                each_misfit_window = misfit_windows[net_sta][category][index]
                each_weight = weights[net_sta][category][index]
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
                total_misfit += (1-each_similarity)*wtotal
                weight_normalize_factor += wtotal
    # normalize
    total_misfit /= weight_normalize_factor
    # get weighted similarity
    weighted_similarity = 1 - total_misfit
    return weighted_similarity


def line_search(alpha_range, t0_range, tau_range, green_virasdf1, green_virasdf2, data_virasdf_body, data_virasdf_surface,
                windows, consider_surface, body_band, surface_band, taper_tmin_tmax, first_arrival_zr, first_arrival_t, baz,
                random_state, tau_raw, t0_raw, stations,  snr_threshold, cc_threshold, deltat_threshold, init_points=9, n_iter=20):
    """
    line_search: given gw, data and essencial information, calculate the optimized parameters combination.
    """
    # do the preprocessing for the virasdf from the green func
    body_green_virasdf1, surface_green_virasdf1 = prepare_green_func(
        green_virasdf1, body_band, surface_band, taper_tmin_tmax, consider_surface)
    body_green_virasdf2, surface_green_virasdf2 = prepare_green_func(
        green_virasdf2, body_band, surface_band, taper_tmin_tmax, consider_surface)
    # get waveform_length and sampling_rate from data
    rep_net_sta = data_virasdf_body.get_waveforms_list()[0]
    rep_st = data_virasdf_body.get_waveforms()[rep_net_sta]["st"]
    rep_tr = rep_st[0]
    waveform_length = rep_tr.stats.endtime - rep_tr.stats.starttime
    sampling_rate = rep_tr.stats.sampling_rate
    # get the weights using tau_raw,t0_raw and alpha=0
    rep_misfit_windows = forward_misfit_windows(tau=tau_raw,
                                                alpha=0,
                                                t0=t0_raw,
                                                body_green_virasdf1=body_green_virasdf1,
                                                body_green_virasdf2=body_green_virasdf2,
                                                body_green_virasdf_raw=body_green_virasdf1,
                                                surface_green_virasdf1=surface_green_virasdf1,
                                                surface_green_virasdf2=surface_green_virasdf2,
                                                surface_green_virasdf_raw=surface_green_virasdf1,
                                                data_virasdf_body=data_virasdf_body,
                                                data_virasdf_surface=data_virasdf_surface,
                                                windows=windows,
                                                waveform_length=waveform_length,
                                                sampling_rate=sampling_rate,
                                                consider_surface=consider_surface,
                                                first_arrival_zr=first_arrival_zr,
                                                first_arrival_t=first_arrival_t,
                                                baz=baz)
    weights = get_weights_for_all(
        rep_misfit_windows, stations, snr_threshold, cc_threshold, deltat_threshold)
    # using the weight to calculate the similarity in the last step:
    weighted_similarity_raw = get_weighted_similarity(
        rep_misfit_windows, weights)
    # wrap forward kernel
    forward = partial(forward_kernel, body_green_virasdf1=body_green_virasdf1,
                      body_green_virasdf2=body_green_virasdf2,
                      body_green_virasdf_raw=body_green_virasdf1,
                      surface_green_virasdf1=surface_green_virasdf1,
                      surface_green_virasdf2=surface_green_virasdf2,
                      surface_green_virasdf_raw=surface_green_virasdf1,
                      data_virasdf_body=data_virasdf_body,
                      data_virasdf_surface=data_virasdf_surface,
                      weights=weights,
                      windows=windows,
                      waveform_length=waveform_length,
                      sampling_rate=sampling_rate,
                      consider_surface=consider_surface,
                      first_arrival_zr=first_arrival_zr,
                      first_arrival_t=first_arrival_t,
                      baz=baz)
    # do the optimization
    pbounds = {'alpha': alpha_range, 'tau': tau_range, 't0': t0_range}
    optimizer = BayesianOptimization(
        f=forward,
        pbounds=pbounds,
        random_state=random_state,
    )
    # probe the central point
    optimizer.probe(
        params={"alpha": 0,
                "tau": tau_raw,
                "t0": t0_raw},
        lazy=True,
    )
    optimizer.maximize(
        init_points=init_points,
        n_iter=n_iter,
    )
    # get the best parameter combination and maximum weighted similarity
    return weighted_similarity_raw, optimizer.max
