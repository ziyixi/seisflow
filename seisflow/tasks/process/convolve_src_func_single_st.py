"""
convolve_src_func_single_st.py: convolve the waveform with the source time function of the green function with the source time function.
"""
import numpy as np
import obspy
import scipy


def source_time_func(tau, dt):
    alpha = 1.628/tau
    # get number of points at each side
    N = int(np.ceil(1.5*tau/dt))
    result = np.zeros(2*N+1)
    # will shift later
    source_t = np.arange(-N, N+1)*dt
    for index, ttime in enumerate(source_t):
        exponentval = alpha**2 * ttime**2
        if(exponentval < 50):
            result[index] = alpha*np.exp(-exponentval)/np.sqrt(np.pi)
        else:
            result[index] = 0
    return result


def conv_sf_and_st(green_st, sf, t0):
    conv_st = obspy.Stream()
    for green_tr in green_st:
        conv_tr = green_tr.copy()
        conv_tr.data = scipy.signal.convolve(
            green_tr.data, sf*green_tr.stats.delta, mode='same', method='auto')
        conv_tr.stats.starttime = conv_tr.stats.starttime+t0
        conv_st += conv_tr
    return conv_st
