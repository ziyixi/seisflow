"""
line_search_source_independent_loc_focal: do line search based on the independent loc/focal information.
"""
from copy import copy
from functools import partial

from bayes_opt import BayesianOptimization

from ..adjoint.calculate_adjoint_source_zerolagcc_one_event import \
    get_weights_for_all
from ..process.convolve_src_func_single_st import source_time_func
from ..windows.calculate_misfit_windows import calculate_misfit_windows
from .line_search_source import (conv_and_postprocess, get_weighted_similarity,
                                 prepare_green_func)


def line_search(alpha_loc_range, alpha_focal_range, t0_range, tau_range, green_virasdf1, green_virasdf2_loc, green_virasdf2_focal, data_virasdf_body, data_virasdf_surface,
                windows, consider_surface, body_band, surface_band, taper_tmin_tmax, first_arrival_zr, first_arrival_t, baz,
                random_state, tau_raw, t0_raw, stations, snr_threshold, cc_threshold, deltat_threshold, init_points=9, n_iter=20):
    # * do preprocessing
    body_green_virasdf1, surface_green_virasdf1 = prepare_green_func(
        green_virasdf1, body_band, surface_band, taper_tmin_tmax, consider_surface)
    body_green_virasdf2_loc, surface_green_virasdf2_loc = prepare_green_func(
        green_virasdf2_loc, body_band, surface_band, taper_tmin_tmax, consider_surface)
    body_green_virasdf2_focal, surface_green_virasdf2_focal = prepare_green_func(
        green_virasdf2_focal, body_band, surface_band, taper_tmin_tmax, consider_surface)
    # * get waveform_length and sampling_rate from data
    rep_net_sta = data_virasdf_body.get_waveforms_list()[0]
    rep_st = data_virasdf_body.get_waveforms()[rep_net_sta]["st"]
    rep_tr = rep_st[0]
    waveform_length = rep_tr.stats.endtime - rep_tr.stats.starttime
    sampling_rate = rep_tr.stats.sampling_rate
    # * get ref misfit_windows (for one predefinied points and get the weight)
    rep_misfit_windows = forward_misfit_windows(tau=tau_raw,
                                                alpha_loc=0,
                                                alpha_focal=0,
                                                t0=t0_raw,
                                                body_green_virasdf1=body_green_virasdf1,
                                                body_green_virasdf2_loc=body_green_virasdf2_loc,
                                                body_green_virasdf2_focal=body_green_virasdf2_focal,
                                                body_green_virasdf_raw=body_green_virasdf1,
                                                surface_green_virasdf1=surface_green_virasdf1,
                                                surface_green_virasdf2_loc=surface_green_virasdf2_loc,
                                                surface_green_virasdf2_focal=surface_green_virasdf2_focal,
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
        rep_misfit_windows, stations, snr_threshold, cc_threshold, deltat_threshold, False)
    # using the weight to calculate the similarity in the last step:
    weighted_similarity_raw = get_weighted_similarity(
        rep_misfit_windows, weights)
    # * wrap forward kernel
    forward = partial(forward_kernel, body_green_virasdf1=body_green_virasdf1,
                      body_green_virasdf2_loc=body_green_virasdf2_loc,
                      body_green_virasdf2_focal=body_green_virasdf2_focal,
                      body_green_virasdf_raw=body_green_virasdf1,
                      surface_green_virasdf1=surface_green_virasdf1,
                      surface_green_virasdf2_loc=surface_green_virasdf2_loc,
                      surface_green_virasdf2_focal=surface_green_virasdf2_focal,
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
    if (tau_range == "fixed"):
        tau_range = (tau_raw, tau_raw)
    pbounds = {'alpha_loc': alpha_loc_range,
               'alpha_focal': alpha_focal_range, 'tau': tau_range, 't0': t0_range}
    optimizer = BayesianOptimization(
        f=forward,
        pbounds=pbounds,
        random_state=random_state,
    )
    # probe the central point
    optimizer.probe(
        params={"alpha_loc": 0,
                "alpha_focal": 0,
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


def get_perturbed_gf(green_virasdf1, green_virasdf2_loc, green_virasdf2_focal, green_virasdf_raw, alpha_loc, alpha_focal):
    """
    get_perturbed_gf: get perturbed asdf
    """
    green_virasdf_perturbed = copy(green_virasdf_raw)
    for each_net_sta in green_virasdf_perturbed.get_waveforms_list():
        st1 = green_virasdf1.get_waveforms()[each_net_sta]["st"]
        st2_loc = green_virasdf2_loc.get_waveforms()[each_net_sta]["st"]
        st2_focal = green_virasdf2_focal.get_waveforms()[each_net_sta]["st"]
        st_perturbed = green_virasdf_perturbed.get_waveforms()[
            each_net_sta]["st"]
        for index in range(len(st_perturbed)):
            st_perturbed[index].data = st_perturbed[index].data + \
                alpha_loc * (st1[index].data - st2_loc[index].data) + \
                alpha_focal * (st1[index].data - st2_focal[index].data)
        green_virasdf_perturbed.update_st(st_perturbed, each_net_sta)
    return green_virasdf_perturbed


def forward_misfit_windows(alpha_loc, alpha_focal, t0, tau, body_green_virasdf1, body_green_virasdf2_loc, body_green_virasdf2_focal, body_green_virasdf_raw, surface_green_virasdf1,
                           surface_green_virasdf2_loc, surface_green_virasdf2_focal, surface_green_virasdf_raw, data_virasdf_body, data_virasdf_surface,
                           windows, waveform_length, sampling_rate, consider_surface, first_arrival_zr, first_arrival_t, baz):
    """
    forward_misfit_windows: get misfit_windows in the forward simulation
    """
    # get perturbed virasdf
    body_green_virasdf_perturbed = get_perturbed_gf(
        body_green_virasdf1, body_green_virasdf2_loc, body_green_virasdf2_focal, body_green_virasdf_raw, alpha_loc, alpha_focal)
    if (consider_surface):
        surface_green_virasdf_perturbed = get_perturbed_gf(
            surface_green_virasdf1, surface_green_virasdf2_loc, surface_green_virasdf2_focal, surface_green_virasdf_raw, alpha_loc, alpha_focal)
    else:
        surface_green_virasdf_perturbed = None
    # get sf
    rep_net_sta = body_green_virasdf_perturbed.get_waveforms_list()[0]
    dt = body_green_virasdf_perturbed.get_waveforms()[
        rep_net_sta]["st"][0].stats.delta
    sf = source_time_func(tau, dt)
    # conv and post process
    body_conv_virasdf_perturbed = conv_and_postprocess(
        body_green_virasdf_perturbed, sf, t0, waveform_length, sampling_rate, True)
    surface_conv_virasdf_perturbed = None
    if (consider_surface):
        surface_conv_virasdf_perturbed = conv_and_postprocess(
            surface_green_virasdf_perturbed, sf, t0, waveform_length, sampling_rate, True)
    # calculate to get misfit windows
    misfit_windows = calculate_misfit_windows(windows, consider_surface, data_virasdf_body, body_conv_virasdf_perturbed,
                                              data_virasdf_surface, surface_conv_virasdf_perturbed, first_arrival_zr, first_arrival_t, baz)
    return misfit_windows


def forward_kernel(alpha_loc, alpha_focal, t0, tau, body_green_virasdf1, body_green_virasdf2_loc, body_green_virasdf2_focal, body_green_virasdf_raw, surface_green_virasdf1,
                   surface_green_virasdf2_loc, surface_green_virasdf2_focal, surface_green_virasdf_raw, data_virasdf_body, data_virasdf_surface,
                   weights, windows, waveform_length, sampling_rate, consider_surface, first_arrival_zr, first_arrival_t, baz):
    """
    forward_kernel: the kernel function for the forward function used in the optimization
    """
    misfit_windows = forward_misfit_windows(alpha_loc, alpha_focal, t0, tau, body_green_virasdf1, body_green_virasdf2_loc, body_green_virasdf2_focal, body_green_virasdf_raw, surface_green_virasdf1,
                                            surface_green_virasdf2_loc, surface_green_virasdf2_focal, surface_green_virasdf_raw, data_virasdf_body, data_virasdf_surface,
                                            windows, waveform_length, sampling_rate, consider_surface, first_arrival_zr, first_arrival_t, baz)
    weighted_similarity = get_weighted_similarity(misfit_windows, weights)
    return weighted_similarity
