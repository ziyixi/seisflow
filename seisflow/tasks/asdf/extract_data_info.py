"""
extract_data_info.py: the core function to extract data information.
"""
from collections import namedtuple

import numpy as np
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees
from obspy.taup import TauPyModel

PHASE_LIST = ["s", "S", "sS", "SS", "p", "P",
              "pP", "sP", "PP", "3.3kmps", "4.6kmps", "ScS"]
EventPair = namedtuple('EventPair', ['gcmtid', 'lat', 'lon', 'dep', 'time'])


def load_station_info(station_fname):
    """
    load_station_info: load stations information.
    """
    # station file
    stations = np.loadtxt(station_fname, dtype=np.str)
    return stations


def kernel(event_pair, stations=None):
    """
    kernel: the core function to get the information.
    """
    model = TauPyModel(model="ak135")
    # for each station, we calculate the travel time
    result = {}
    evla = event_pair.lat
    evlo = event_pair.lon
    evdp = event_pair.dep/1000
    event_time = event_pair.time

    for row in stations:
        result_template = {
            "event_time": event_time,
            "evla": evla,
            "evlo": evlo,
            "evdp": evdp,
            "gcarc": None,
            "az": None,
            "baz": None,
            "S": None,
            "sS": None,
            "SS": None,
            "P": None,
            "pP": None,
            "sP": None,
            "PP": None,
            "3.3kmps": None,
            "4.6kmps": None,
            "ScS": None
        }
        station = row[0]
        network = row[1]
        net_sta = f"{network}.{station}"
        stla = float(row[2])
        stlo = float(row[3])
        arrivals = model.get_travel_times_geo(
            evdp, evla, evlo, stla, stlo, PHASE_LIST)

        gcarc_m, az, baz = gps2dist_azimuth(evla, evlo, stla, stlo)
        gcarc = kilometers2degrees(gcarc_m / 1000)
        result_template["gcarc"] = gcarc
        result_template["az"] = az
        result_template["baz"] = baz

        for each_arrival in arrivals:
            name = each_arrival.name
            time = each_arrival.time
            if ((name in PHASE_LIST)):
                if (name == "p"):
                    name = "P"
                if (name == "s"):
                    name = "S"
                if (result_template[name] == None):
                    result_template[name] = time
        result[net_sta] = result_template
    return result


def extract_data_info(gcmtid, evla, evlo, evdp, event_time, stations):
    used_event_pair = EventPair(
        gcmtid=gcmtid, lat=evla, lon=evlo, dep=evdp, time=event_time)
    data_info = kernel(used_event_pair, stations=stations)
    return data_info
