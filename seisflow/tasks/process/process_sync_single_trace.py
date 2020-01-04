"""
process_sync_single_trace.py: process only a single trace of sync, used by other scripts.
"""

import numpy as np
import obspy
from obspy.signal.util import _npts2nfft


def sync_remove_response(pre_filt, st):
    """
    mimic obspy.remove_response, but only do the frequency taper
    """
    obspy.core.util.misc.limit_numpy_fft_cache()
    for trace in st:
        data = trace.data.astype(np.float64)
        npts = len(data)
        nfft = _npts2nfft(npts)
        data = np.fft.rfft(data, n=nfft)
        t_samp = trace.stats.delta
        fy = 1 / (t_samp * 2.0)
        freqs = np.linspace(0, fy, nfft // 2 + 1).astype(np.float64)
        freq_domain_taper = obspy.signal.invsim.cosine_sac_taper(
            freqs, flimit=pre_filt)
        data *= freq_domain_taper
        data = np.fft.irfft(data)[0:npts]
        trace.data = data

    return st


def check_time(st, event_time, waveform_length):
    for trace in st:
        starttime = trace.stats.starttime
        endtime = trace.stats.endtime
        if(starttime > event_time):
            return -1
        if(endtime < event_time+waveform_length):
            return -1
    return 0


def process_sync_single_trace(st, event_time, waveform_length, taper_tmin_tmax, min_period, max_period, sampling_rate):
    status_code = check_time(st, event_time, waveform_length)
    if(status_code == 0):
        pass
    elif(status_code == -1):
        return
    else:
        raise Exception("unknown status code")
    st.trim(event_time, event_time+waveform_length)
    st.detrend("demean")
    st.detrend("linear")
    st.taper(max_percentage=0.05, type="hann")

    tmin, tmax = map(float, taper_tmin_tmax.split(","))
    f2 = 1.0 / tmax
    f3 = 1.0 / tmin
    f1 = 0.5 * f2
    f4 = 2.0 * f3
    pre_filt = (f1, f2, f3, f4)
    sync_remove_response(pre_filt, st)

    st.interpolate(sampling_rate=sampling_rate)

    # bandpass filter
    st.filter("bandpass", freqmin=1.0/max_period,
              freqmax=1.0/min_period, corners=2, zerophase=True)
    return st
