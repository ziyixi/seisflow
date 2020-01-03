"""
tao_2018_ggg_windows.py: get windows described in tao et al. ggg, 2018. Only works for one event (gcmtid).
"""

import obspy

from .window import Window, Windows_collection

phases = ["S", "sS", "SS", "P", "pP",
          "sP", "PP", "3.3kmps", "4.6kmps", "ScS"]


def generate_windows(gcmtid, traveltime, event_time, time_length):
    # get all the net_sta
    all_net_sta = list(traveltime.keys())
    # for each net_sta, get all the windows, and merge them accordingly.
    phases_zr = ["P", "pP", "sP", "PP", "S", "sS", "SS"]
    phases_t = ["ScS", "S", "sS", "SS"]
    result = {}
    for each_net_sta in all_net_sta:
        result[each_net_sta] = {
            "z": Windows_collection(),
            "r": Windows_collection(),
            "t": Windows_collection(),
            "surface_z": Windows_collection(),
            "surface_r": Windows_collection(),
            "surface_t": Windows_collection()
        }
        # z
        for each_phase in phases_zr:
            each_phase_traveltime = traveltime[each_net_sta][each_phase]
            if (each_phase_traveltime == None):
                continue
            elif (each_phase_traveltime > time_length):
                continue
            else:
                left = event_time + each_phase_traveltime - 20
                if (left < event_time):
                    left = event_time
                right = event_time + each_phase_traveltime + 50
                if (right > event_time + time_length):
                    right = event_time + time_length
                channel = "Z"
                network = each_net_sta.split(".")[0]
                station = each_net_sta.split(".")[1]
                phases = [each_phase]
                result[each_net_sta]["z"].append_window(Window(
                    left=left, right=right, channel=channel, network=network, station=station, phases=phases, gcmtid=gcmtid))
        # r
        for each_phase in phases_zr:
            each_phase_traveltime = traveltime[each_net_sta][each_phase]
            if (each_phase_traveltime == None):
                continue
            elif (each_phase_traveltime > time_length):
                continue
            else:
                left = event_time + each_phase_traveltime - 20
                if (left < event_time):
                    left = event_time
                right = event_time + each_phase_traveltime + 50
                if (right > event_time + time_length):
                    right = event_time + time_length
                channel = "R"
                network = each_net_sta.split(".")[0]
                station = each_net_sta.split(".")[1]
                phases = [each_phase]
                result[each_net_sta]["r"].append_window(Window(
                    left=left, right=right, channel=channel, network=network, station=station, phases=phases, gcmtid=gcmtid))
        # t
        for each_phase in phases_t:
            each_phase_traveltime = traveltime[each_net_sta][each_phase]
            if (each_phase_traveltime == None):
                continue
            elif (each_phase_traveltime > time_length):
                continue
            else:
                left = event_time + each_phase_traveltime - 20
                if (left < event_time):
                    left = event_time
                right = event_time + each_phase_traveltime + 50
                if (right > event_time + time_length):
                    right = event_time + time_length
                channel = "T"
                network = each_net_sta.split(".")[0]
                station = each_net_sta.split(".")[1]
                phases = [each_phase]
                result[each_net_sta]["t"].append_window(Window(
                    left=left, right=right, channel=channel, network=network, station=station, phases=phases, gcmtid=gcmtid))
        # surface_z
        phase_travel_time_left = traveltime[each_net_sta]["4.6kmps"]
        phase_travel_time_right = traveltime[each_net_sta]["3.3kmps"]
        if(phase_travel_time_right < time_length):
            left = event_time+phase_travel_time_left-50
            if (left < event_time):
                left = event_time
            right = event_time + phase_travel_time_right + 50
            if (right > event_time + time_length):
                right = event_time+time_length
            channel = "Z"
            network = each_net_sta.split(".")[0]
            station = each_net_sta.split(".")[1]
            phases = ["surface_z"]
            result[each_net_sta]["surface_z"].append_window(Window(
                left=left, right=right, channel=channel, network=network, station=station, phases=phases, gcmtid=gcmtid))
        # surface_r
        phase_travel_time_left = traveltime[each_net_sta]["4.6kmps"]
        phase_travel_time_right = traveltime[each_net_sta]["3.3kmps"]
        if(phase_travel_time_right < time_length):
            left = event_time+phase_travel_time_left-50
            if (left < event_time):
                left = event_time
            right = event_time + phase_travel_time_right + 50
            if (right > event_time + time_length):
                right = event_time+time_length
            channel = "R"
            network = each_net_sta.split(".")[0]
            station = each_net_sta.split(".")[1]
            phases = ["surface_r"]
            result[each_net_sta]["surface_r"].append_window(Window(
                left=left, right=right, channel=channel, network=network, station=station, phases=phases, gcmtid=gcmtid))
        # surface_t
        phase_travel_time_left = traveltime[each_net_sta]["4.6kmps"]
        phase_travel_time_right = traveltime[each_net_sta]["3.3kmps"]
        if(phase_travel_time_right < time_length):
            left = event_time+phase_travel_time_left-50
            if (left < event_time):
                left = event_time
            right = event_time + phase_travel_time_right + 50
            if (right > event_time + time_length):
                right = event_time+time_length
            channel = "T"
            network = each_net_sta.split(".")[0]
            station = each_net_sta.split(".")[1]
            phases = ["surface_t"]
            result[each_net_sta]["surface_t"].append_window(Window(
                left=left, right=right, channel=channel, network=network, station=station, phases=phases, gcmtid=gcmtid))

    # at the moment, we have the result dict: gcmtid->net_sta->type
    # merge the windows
        for each_net_sta in result:
            result_each_net_sta = result[each_net_sta]
            result_each_net_sta["z"].merge_windows()
            result_each_net_sta["r"].merge_windows()
            result_each_net_sta["t"].merge_windows()

    return result
