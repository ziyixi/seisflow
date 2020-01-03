"""
process sync file. 
"""
from os.path import join

# fix a bug in intel
import mpi4py
import numpy as np
import obspy
import pyasdf
from mpi4py import MPI
from obspy.signal.util import _npts2nfft

mpi4py.rc.recv_mprobe = False


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
isroot = (rank == 0)


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


def check_time(st, event_time, waveform_length, inv):
    for trace in st:
        starttime = trace.stats.starttime
        endtime = trace.stats.endtime
        if(starttime > event_time):
            return -1
        if(endtime < event_time+waveform_length):
            return -1
    return 0


def process_single_event(min_periods, max_periods, taper_tmin_tmax, asdf_filename, waveform_length, sampling_rate, output_directory):
    tmin, tmax = map(float, taper_tmin_tmax.split(","))
    # with pyasdf.ASDFDataSet(asdf_filename) as ds:
    with pyasdf.ASDFDataSet(asdf_filename, mode="r") as ds:
        # some parameters
        event = ds.events[0]
        origin = event.preferred_origin() or event.origins[0]
        event_time = origin.time
        event_latitude = origin.latitude
        event_longitude = origin.longitude

        for min_period, max_period in zip(min_periods, max_periods):
            f2 = 1.0 / tmax
            f3 = 1.0 / tmin
            f1 = 0.5 * f2
            f4 = 2.0 * f3
            pre_filt = (f1, f2, f3, f4)
            # inv will be None if no inv provided

            def process_function(st, inv):
                status_code = check_time(st, event_time, waveform_length, inv)
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

                sync_remove_response(pre_filt, st)

                st.interpolate(sampling_rate=sampling_rate)

                # bandpass filter
                st.filter("bandpass", freqmin=1.0/max_period,
                          freqmax=1.0/min_period, corners=2, zerophase=True)

                # Convert to single precision to save space.
                for tr in st:
                    tr.data = np.require(tr.data, dtype="float32")

                return st
            tag_name = "preprocessed_%is_to_%is" % (
                int(min_period), int(max_period))
            tag_map = {
                "synthetic": tag_name
            }

            output_name_head = asdf_filename.split("/")[-1].split(".")[0]
            ds.process(process_function, join(
                output_directory, output_name_head+"."+tag_name + ".h5"), tag_map=tag_map)
