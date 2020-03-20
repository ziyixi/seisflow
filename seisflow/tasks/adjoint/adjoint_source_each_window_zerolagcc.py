import numpy as np


def calculate_adjoint_source_each_window(mistit_window, raw_sync_asdf_trace, sync_asdf_trace, data_asdf_trace, mintime, maxtime, raw_sync_asdf_trace_adjust_time=0):
    # firstly we prepare the windowed trace with tapering.
    # all the input traces should be a copy
    # since we are not going to invert for the half duration, using raw_sync_asdf_trace calculated using the raw CMTSOLUTION is appropriate for the reference purpose.
    similarity = mistit_window.similarity

    # hope they will have the same length, or we have to do some trick
    # actually they don't have the same length, we might do some trick here
    windowed_data_trace = data_asdf_trace.slice(
        mistit_window.left, mistit_window.right)
    windowed_sync_trace = sync_asdf_trace.slice(
        mistit_window.left, mistit_window.right)

    true_windowed_sync_time = windowed_sync_trace.stats.starttime
    offset_window_unwindowed = round(
        (true_windowed_sync_time - sync_asdf_trace.stats.starttime) / sync_asdf_trace.stats.delta)
    # ! fix a bug here, note there is possibility that the event time of raw_sync_asdf_trace and sync_asdf_trace is different,
    # ! which means sync_asdf_trace is not directly from raw_sync_asdf_trace

    # ! note here we introduce another bug, if offset_unwindowed_raw<0 due to raw_sync_asdf_trace_adjust_time, there will be a problem
    # ! to fix it,
    offset_unwindowed_raw = round(
        (sync_asdf_trace.stats.starttime-(raw_sync_asdf_trace.stats.starttime+raw_sync_asdf_trace_adjust_time))/(raw_sync_asdf_trace.stats.delta))

    # to keep the windowed data and the windowed sync have the same npts
    windowed_sync_trace.stats.starttime = windowed_data_trace.stats.starttime
    newright = min(windowed_sync_trace.stats.endtime,
                   windowed_data_trace.stats.endtime)
    windowed_data_trace = windowed_data_trace.slice(windowed_data_trace.stats.starttime, newright
                                                    )
    windowed_sync_trace = windowed_sync_trace.slice(windowed_data_trace.stats.starttime, newright
                                                    )

    # taper
    windowed_data_trace.taper(0.05, type="hann")
    windowed_sync_trace.taper(0.05, type="hann")
    obs_norm = np.sqrt(np.sum(windowed_data_trace.data ** 2))
    syn_norm = np.sqrt(np.sum(windowed_sync_trace.data ** 2))
    Nw = obs_norm * syn_norm
    Aw = similarity * obs_norm / syn_norm
    syn_delta = windowed_sync_trace.stats.delta
    raw_syn_sampling_rate = raw_sync_asdf_trace.stats.sampling_rate

    # here the adjoint_source_windowed is the windowed trace
    adjoint_source_windowed = windowed_sync_trace.copy()
    adjoint_source_windowed.data[:] = 0.0
    obs_filt_win = windowed_data_trace.data
    syn_filt_win = windowed_sync_trace.data
    # !
    # ! note, here is an error. In tao's code, he is using cc as the so called misfit. Thus comparing with the paper,
    # ! there will be no negative for the below formula. I am using the misfit, thus should add a minus. Fixed here!!
    # !
    adjoint_source_windowed.data = -(
        obs_filt_win - Aw * syn_filt_win) / Nw / syn_delta
    # apply taper and bandpass
    adjoint_source_windowed.taper(0.05, type="hann")
    adjoint_source_windowed.filter(
        "bandpass", freqmin=1/maxtime, freqmax=1/mintime, corners=2, zerophase=True)
    # retrive to unwindowed trace
    len_window = len(adjoint_source_windowed.data)
    adjoint_source_unwindowed = sync_asdf_trace.copy()
    adjoint_source_unwindowed.data[:] = 0.0
    adjoint_source_unwindowed.data[offset_window_unwindowed:
                                   offset_window_unwindowed+len_window] = adjoint_source_windowed.data
    # interpolate to the raw sync
    adjoint_source_unwindowed.interpolate(
        sampling_rate=raw_syn_sampling_rate)
    len_adjoint_source_unwindowed = len(adjoint_source_unwindowed.data)
    adjoint_source_final = raw_sync_asdf_trace.copy()
    adjoint_source_final.data[:] = 0.0
    # * fix rhe bug when offset_unwindowed_raw<0
    if(offset_unwindowed_raw >= 0):
        adjoint_source_final.data[offset_unwindowed_raw:
                                  offset_unwindowed_raw + len_adjoint_source_unwindowed] = adjoint_source_unwindowed.data
    else:
        # we assume [-offset,0] has no adjoint source
        offset = -offset_unwindowed_raw
        adjoint_source_final.data[0:
                                  offset_unwindowed_raw + len_adjoint_source_unwindowed] = adjoint_source_unwindowed.data[offset:]
    # * fix the problem of nan, in case the data value are all 0.
    if (np.isnan(mistit_window.snr_energy)):
        adjoint_source_final.data[:] = 0.0
    return adjoint_source_final
