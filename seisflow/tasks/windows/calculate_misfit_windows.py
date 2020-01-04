"""
calculate_misfit_windows.py: get the misfit window given the asdf files and the windows.
"""
import obspy
import pyasdf

from .misfit_window import Misfit_window
from .window import Window, Windows_collection


def get_used_net_sta(windows, data_asdf_body_path):
    net_sta_sync = list(windows.keys())
    data_asdf_body = pyasdf.ASDFDataSet(data_asdf_body_path, mode="r")
    net_sta_data = data_asdf_body.waveforms.list()
    used_net_sta = list(set(net_sta_data) & set(net_sta_sync))
    del data_asdf_body
    return used_net_sta


def prepare_windows(windows_used_event,  consider_surface, used_net_sta):
    """
    Generate misfit windows and split different component from the windows.
    """
    new_windows = {}
    for net_sta in used_net_sta:
        if(consider_surface):
            new_windows[net_sta] = {
                "z": Windows_collection(),
                "r": Windows_collection(),
                "t": Windows_collection(),
                "surface_z": Windows_collection(),
                "surface_r": Windows_collection(),
                "surface_t": Windows_collection()
            }
        else:
            new_windows[net_sta] = {
                "z": Windows_collection(),
                "r": Windows_collection(),
                "t": Windows_collection()
            }
        old_windows_z = windows_used_event[net_sta]["z"].windows
        old_windows_r = windows_used_event[net_sta]["r"].windows
        old_windows_t = windows_used_event[net_sta]["t"].windows
        if(consider_surface):
            old_windows_surface_z = windows_used_event[net_sta]["surface_z"].windows
            old_windows_surface_r = windows_used_event[net_sta]["surface_r"].windows
            old_windows_surface_t = windows_used_event[net_sta]["surface_t"].windows
        # update to the new_windows
        for each_window in old_windows_z:
            new_windows[net_sta]["z"].append_window(
                Misfit_window(each_window))
        for each_window in old_windows_r:
            new_windows[net_sta]["r"].append_window(
                Misfit_window(each_window))
        for each_window in old_windows_t:
            new_windows[net_sta]["t"].append_window(
                Misfit_window(each_window))
        if (consider_surface):
            for each_window in old_windows_surface_z:
                new_windows[net_sta]["surface_z"].append_window(
                    Misfit_window(each_window))
            for each_window in old_windows_surface_r:
                new_windows[net_sta]["surface_r"].append_window(
                    Misfit_window(each_window))
            for each_window in old_windows_surface_t:
                new_windows[net_sta]["surface_t"].append_window(
                    Misfit_window(each_window))
    return new_windows


def calculate_snr_cc_deltat(data_asdf_body_path, sync_asdf_body_path, data_asdf_surface_path, sync_asdf_surface_path, misfit_windows, first_arrival_zr, first_arrival_t, baz, consider_surface):
    data_asdf_body = pyasdf.ASDFDataSet(data_asdf_body_path)
    sync_asdf_body = pyasdf.ASDFDataSet(sync_asdf_body_path)
    if(consider_surface):
        data_asdf_surface = pyasdf.ASDFDataSet(data_asdf_surface_path)
        sync_asdf_surface = pyasdf.ASDFDataSet(sync_asdf_surface_path)
    for net_sta in misfit_windows:
        for category in misfit_windows[net_sta]:
            for each_window in misfit_windows[net_sta][category].windows:
                # update the first arrival
                if(each_window.channel == "Z" or each_window.channel == "R"):
                    each_window.update_first_arrival_baz(
                        first_arrival_zr, baz)
                elif (each_window.channel == "T"):
                    each_window.update_first_arrival_baz(
                        first_arrival_t, baz)
                else:
                    raise Exception(
                        "channel is not correct in updating the first arrival!")
                # update snr,deltat and cc
                if((category == "z") or (category == "r") or (category == "t")):
                    each_window.update_snr(data_asdf_body, sync_asdf_body)
                    each_window.update_cc_deltat(
                        data_asdf_body, sync_asdf_body)
                elif((category == "surface_z") or (category == "surface_r") or (category == "surface_t")):
                    if (consider_surface):
                        each_window.update_snr(
                            data_asdf_surface, sync_asdf_body)
                        each_window.update_cc_deltat(
                            data_asdf_surface, sync_asdf_surface)
                else:
                    raise Exception(
                        "category is not correct in calculatng snr,delta and cc")
    del data_asdf_body
    del sync_asdf_body
    if (consider_surface):
        del data_asdf_surface
        del sync_asdf_surface
    return misfit_windows


def calculate_misfit_windows(windows, consider_surface, data_asdf_body_path, sync_asdf_body_path, data_asdf_surface_path, sync_asdf_surface_path, first_arrival_zr, first_arrival_t, baz):
    """
    calculate_misfit_windows: calculate misfit windows.
    """
    used_net_sta = get_used_net_sta(windows, data_asdf_body_path)
    misfit_windows = prepare_windows(windows, consider_surface, used_net_sta)
    misfit_windows = calculate_snr_cc_deltat(data_asdf_body_path, sync_asdf_body_path, data_asdf_surface_path,
                                             sync_asdf_surface_path, misfit_windows, first_arrival_zr, first_arrival_t, baz, consider_surface)
    return misfit_windows
