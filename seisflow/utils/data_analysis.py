"""
data_analysis.py: tools to analyze data.
"""
import numpy as np


def get_windows_deltat(misfit_windows, exclude_p, exclude_s, component, snr_threshold, phase, all_threshold=False, cc_threshold=None, deltat_threshold=None):
    result = []
    rep_net_sta = sorted(misfit_windows.keys())[0]
    if (component not in misfit_windows[rep_net_sta]):
        return np.array([])
    for net_sta in sorted(misfit_windows.keys()):
        for each_window in misfit_windows[net_sta][component].windows:
            condition = True
            value = each_window.deltat
            if(exclude_p and "P" in each_window.phases):
                condition = False
            if (exclude_s and "S" in each_window.phases):
                condition = False
            if (phase not in each_window.phases):
                condition = False
            if (condition):
                if ((not np.isnan(each_window.snr_energy)) and (not np.isnan(each_window.snr_amp)) and (not np.isnan(each_window.deltat))
                        and (not np.isnan(each_window.cc)) and (not np.isnan(each_window.similarity))):
                    if(not all_threshold):
                        if(each_window.snr_energy >= snr_threshold):
                            result.append(value)
                    else:
                        if ((each_window.snr_energy >= snr_threshold) and (each_window.cc >= cc_threshold) and (np.abs(each_window.deltat) <= deltat_threshold)):
                            result.append(value)
    return np.array(result)


def get_windows_snr_energy(misfit_windows, exclude_p, exclude_s, component, snr_threshold, phase):
    result = []
    for net_sta in sorted(misfit_windows.keys()):
        for each_window in misfit_windows[net_sta][component].windows:
            condition = True
            value = each_window.snr_energy
            if(exclude_p and "P" in each_window.phases):
                condition = False
            if (exclude_s and "S" in each_window.phases):
                condition = False
            if (phase not in each_window.phases):
                condition = False
            if (condition):
                if ((not np.isnan(each_window.snr_energy)) and (not np.isnan(each_window.snr_amp)) and (not np.isnan(each_window.deltat))
                        and (not np.isnan(each_window.cc)) and (not np.isnan(each_window.similarity))):
                    if(each_window.snr_energy > snr_threshold):
                        result.append(value)
    return np.array(result)


def get_windows_cc(misfit_windows, exclude_p, exclude_s, component, snr_threshold, phase, all_threshold=False, cc_threshold=None, deltat_threshold=None):
    result = []
    rep_net_sta = sorted(misfit_windows.keys())[0]
    if (component not in misfit_windows[rep_net_sta]):
        return np.array([])
    for net_sta in sorted(misfit_windows.keys()):
        for each_window in misfit_windows[net_sta][component].windows:
            condition = True
            value = each_window.cc
            if(exclude_p and "P" in each_window.phases):
                condition = False
            if (exclude_s and "S" in each_window.phases):
                condition = False
            if (phase not in each_window.phases):
                condition = False
            if (condition):
                if ((not np.isnan(each_window.snr_energy)) and (not np.isnan(each_window.snr_amp)) and (not np.isnan(each_window.deltat))
                        and (not np.isnan(each_window.cc)) and (not np.isnan(each_window.similarity))):
                    if(not all_threshold):
                        if(each_window.snr_energy >= snr_threshold):
                            result.append(value)
                    else:
                        if ((each_window.snr_energy >= snr_threshold) and (each_window.cc >= cc_threshold) and (np.abs(each_window.deltat) <= deltat_threshold)):
                            result.append(value)
    return np.array(result)


def get_windows_similarity(misfit_windows, exclude_p, exclude_s, component, snr_threshold, phase, all_threshold=False, cc_threshold=None, deltat_threshold=None):
    result = []
    rep_net_sta = sorted(misfit_windows.keys())[0]
    if (component not in misfit_windows[rep_net_sta]):
        return np.array([])
    for net_sta in sorted(misfit_windows.keys()):
        for each_window in misfit_windows[net_sta][component].windows:
            condition = True
            value = each_window.similarity
            if(exclude_p and "P" in each_window.phases):
                condition = False
            if (exclude_s and "S" in each_window.phases):
                condition = False
            if (phase not in each_window.phases):
                condition = False
            if (condition):
                if ((not np.isnan(each_window.snr_energy)) and (not np.isnan(each_window.snr_amp)) and (not np.isnan(each_window.deltat))
                        and (not np.isnan(each_window.cc)) and (not np.isnan(each_window.similarity))):
                    if(not all_threshold):
                        if(each_window.snr_energy >= snr_threshold):
                            result.append(value)
                    else:
                        if ((each_window.snr_energy >= snr_threshold) and (each_window.cc >= cc_threshold) and (np.abs(each_window.deltat) <= deltat_threshold)):
                            result.append(value)
    return np.array(result)


def get_windows_net_sta(misfit_windows, exclude_p, exclude_s, component, snr_threshold, phase):
    result = []
    for net_sta in sorted(misfit_windows.keys()):
        for each_window in misfit_windows[net_sta][component].windows:
            condition = True
            value = net_sta
            if(exclude_p and "P" in each_window.phases):
                condition = False
            if (exclude_s and "S" in each_window.phases):
                condition = False
            if (phase not in each_window.phases):
                condition = False
            if (condition):
                if ((not np.isnan(each_window.snr_energy)) and (not np.isnan(each_window.snr_amp)) and (not np.isnan(each_window.deltat))
                        and (not np.isnan(each_window.cc)) and (not np.isnan(each_window.similarity))):
                    if(each_window.snr_energy > snr_threshold):
                        result.append(value)
    return np.array(result)
