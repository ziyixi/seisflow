"""
convert_green2processed.py: do the convolution and process the green sync.
"""
from .process_sync_single_st import ahead_process_green_function_st, post_process_green_function_st
from .convolve_src_func_single_st import source_time_func, conv_sf_and_st
import copy


def conv_process_single_virasdf(input_virasdf, waveform_length, tau, t0,
                                taper_tmin_tmax, min_period, max_period, sampling_rate):
    """
    convert the green function virasdf to the normal form of asdf.
    """
    # firstly we generate the sf
    rep_net_sta = input_virasdf.get_waveforms_list()[0]
    rep_tr = input_virasdf.get_waveforms()[rep_net_sta]["st"][0]
    dt = rep_tr.stats.delta
    sf = source_time_func(tau, dt)
    event = input_virasdf.events[0]
    origin = event.preferred_origin() or event.origins[0]
    event_time = origin.time
    output_virasdf = copy.copy(input_virasdf)
    for each_net_sta in output_virasdf.get_waveforms_list():
        green_st = output_virasdf.get_waveforms()[each_net_sta]
        green_st = ahead_process_green_function_st(
            green_st, taper_tmin_tmax, min_period, max_period)
        conv_st = conv_sf_and_st(green_st, sf, t0)
        conv_st = post_process_green_function_st(
            conv_st, event_time, waveform_length, sampling_rate)
        output_virasdf.update_st(conv_st, each_net_sta)
    return output_virasdf
