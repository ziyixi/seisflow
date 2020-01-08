"""
process asdf data file.
"""

from os.path import join

# fix a bug in intel
import mpi4py
import numpy as np
import obspy
import pandas as pd
import pyasdf
from mpi4py import MPI
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.invsim import simulate_seismometer
from pyasdf import ASDFDataSet

mpi4py.rc.recv_mprobe = False


# global parameters
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
isroot = (rank == 0)

# CEA_NETWORKS
CEA_NETWORKS = ["AH", "BJ", "BU", "CQ", "FJ", "GD", "GS", "GX", "GZ", "HA", "HB", "HE", "HI", "HL", "HN",
                "JL", "JS", "JX", "LN", "NM", "NX", "QH", "SC", "SD", "SH", "SN", "SX", "TJ", "XJ", "XZ", "YN", "ZJ"]


def check_st_numberlap(st, inv):
    """
    detect overlapping
    """
    if(len(st) == 0):
        return -1
    elif(len(st) < 3):
        return -1
    elif(len(st) == 3):
        channel_set = set()
        for item in st:
            channel_set.add(item.id[-1])
        if(len(channel_set) == 3):
            return 0
        else:
            return -1
    else:
        channel_set = set()
        for item in st:
            channel_set.add(item.id[-1])
        if(len(channel_set) == 3):
            return 1
        elif(len(channel_set) < 3):
            return -1
        else:
            return -1


def modify_time(time_str):
    if(type(time_str) != str):
        return obspy.UTCDateTime("2099-09-01")
    else:
        return obspy.UTCDateTime(time_str)


def func_correct_cea(baz, inv, event_time, correction_data):
    network, station = inv.get_contents()["channels"][0].split(".")[:2]

    # after 2013-09-01T11:52:00Z, all the stations have been corrected
    trunc_datetime = obspy.UTCDateTime("2013-09-01T11:52:00Z")
    if(event_time > trunc_datetime):
        return baz
    else:
        info_for_this_station = correction_data.loc[(
            correction_data.network == network) & (correction_data.station == station) & (correction_data.starttime <= event_time) & (correction_data.endtime >= event_time)]
        if(len(info_for_this_station) == 0):
            return None
        elif(len(info_for_this_station) == 1):
            median_value = info_for_this_station["median"].values[0]
            if(np.isnan(median_value)):
                if_has_been_corrected = (
                    info_for_this_station["endtime"].values[0] == obspy.UTCDateTime("2099-09-01"))
                if(if_has_been_corrected):
                    return baz
                else:
                    return None
            return baz-median_value
        else:
            return None


def check_time(st, event_time, waveform_length, inv):
    for trace in st:
        starttime = trace.stats.starttime
        endtime = trace.stats.endtime
        # add some tolerance here (1 min)
        if(starttime-60 > event_time):
            return -1
        if(endtime < event_time+waveform_length):
            return -1
    return 0


def filter_st(st, inv):
    # we should assure only to add one type of data, as the order of HH,BH,SH (we don't consider the case like
    # only 2 HH but 3 BH.)
    newst = obspy.Stream()
    # get band code status
    band_code = None
    band_code_list = []
    for trace in st:
        theid = trace.id
        net, sta, loc, cha = theid.split(".")
        band_code_list.append(cha[:2])
    if(len(band_code_list) == 0):
        return newst
    else:
        if("HH" in band_code_list):
            band_code = "HH"
        elif("BH" in band_code_list):
            band_code = "BH"
        elif("SH" in band_code_list):
            band_code = "SH"
        else:
            return newst
    for trace in st:
        theid = trace.id
        net, sta, loc, cha = theid.split(".")

        con1 = ((loc == "") or (loc == "00"))
        con2 = (cha[:2] == band_code)
        if(con1 and con2):
            newst += trace
    return newst


def remove_response_paz(st, paz_path, pre_filt):
    """
    remove response using paz file
    """
    for trace in st:
        # get key
        network = trace.stats.network
        station = trace.stats.station
        channel = trace.stats.channel
        key = f"{network}.{station}.{channel}"
        paz_path = join(paz_path, key)
        # try to read paz
        try:
            obspy.io.sac.sacpz.attach_paz(trace, paz_path)
        except:
            return None
        # ndarray
        data = simulate_seismometer(
            trace.data, trace.stats.sampling_rate, paz_remove=trace.stats.paz, water_level=6e9, zero_mean=False, taper=False, pre_filt=pre_filt, sacsim=True)
        trace.data = data
    return st


def process_single_event(min_periods, max_periods, taper_tmin_tmax, asdf_filename, waveform_length, sampling_rate, output_directory, correct_cea, cea_correction_file, paz_path):
    tmin, tmax = map(float, taper_tmin_tmax.split(","))
    # with pyasdf.ASDFDataSet(asdf_filename) as ds:
    ds = pyasdf.ASDFDataSet(asdf_filename, mode="r")

    # load cea correction file
    if(correct_cea):
        correction_data = pd.read_csv(cea_correction_file, sep="|", comment="#", names=[
            "network", "station", "eventno", "mean", "std", "median", "mad", "starttime", "endtime"])
        correction_data["starttime"] = correction_data["starttime"].apply(
            modify_time)
        correction_data["endtime"] = correction_data["endtime"].apply(
            modify_time)

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

        def process_function(st, inv):
            # there are possibility that some stations has multiple loc codes or use HH stations. (should avoid in the future)
            st = filter_st(st, inv)

            # overlap the previous trace
            status_code = check_st_numberlap(st, inv)
            if(status_code == -1):
                return
            elif(status_code == 0):
                pass
            elif(status_code == 1):
                # merge may have roblem (samplign rate is not equal)
                try:
                    st.merge(method=1, fill_value=0, interpolation_samples=0)
                except:
                    return
            else:
                raise Exception("unknown status code")

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

            # st.remove_response(output="DISP", pre_filt=pre_filt, zero_mean=False,
            #                    taper=False, inventory=inv, water_level=None)
            st = remove_response_paz(st, paz_path, pre_filt)
            # the same of removing response with sac
            st.detrend("demean")
            st.detrend("linear")

            if(st == None):
                return
            # sac has already removed mean after removing the response.

            st.interpolate(sampling_rate=sampling_rate)

            station_latitude = inv[0][0].latitude
            station_longitude = inv[0][0].longitude

            # baz is calculated using station and event's location
            # for cea stations, we can directly add an angle to it
            _, baz, _ = gps2dist_azimuth(station_latitude, station_longitude,
                                         event_latitude, event_longitude)

            network = inv.get_contents()['networks'][0]
            if(correct_cea and (network in CEA_NETWORKS)):
                baz = func_correct_cea(
                    baz, inv, event_time, correction_data)
            if(baz == None):
                return

            # we have to limit baz to be in [0,360)
            baz = np.mod(baz, 360)

            components = [tr.stats.channel[-1] for tr in st]
            if "N" in components and "E" in components:
                # there may be some problem in rotating (time span is not equal for three channels)
                try:
                    st.rotate(method="NE->RT", back_azimuth=baz)
                except:
                    return
            else:
                return
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
            "raw": tag_name
        }
        output_name_head = asdf_filename.split("/")[-1].split(".")[0]
        ds.process(process_function, join(
            output_directory, output_name_head+"."+tag_name + ".h5"), tag_map=tag_map)

    del ds
