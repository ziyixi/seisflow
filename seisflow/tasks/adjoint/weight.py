"""
weight.py: calculate weight for each windowed trace.
"""
import numpy as np
from obspy.geodetics import locations2degrees


def cal_cos_weight(value, value1, value2):
    result = 0
    if(value < value1):
        result = 0
    elif(value1 <= value < value2):
        result = 0.5+0.5*np.cos(np.pi*(value2-value)/(value2-value1))
    else:
        result = 1
    if(np.isnan(value)):
        result = 0
    return result


def cal_cos_weight_ops(value, value1, value2):
    result = 0
    if (value < value1):
        result = 1
    elif (value1 <= value < value2):
        result = 0.5 - 0.5 * \
            np.cos(np.pi * (value2 - value) / (value2 - value1))
    else:
        result = 0
    if(np.isnan(value)):
        result = 0
    return result

##################################################################################
# weights used for single windowed trace
##################################################################################


def cal_snr_weight(used_misfit_window, value1, value2):
    """
    cal_snr_weight: use snr_energy to get a weight.
    """
    value = used_misfit_window.snr_energy
    weight = cal_cos_weight(value, value1, value2)
    return weight


def cal_cc_weight(used_misfit_window, value1, value2):
    """
    cal_cc_weight: use cc to get a weight.
    """
    value = used_misfit_window.cc
    weight = cal_cos_weight(value, value1, value2)
    return weight


def cal_deltat_weight(used_misfit_window, value1, value2):
    """
    cal_deltat_weight: use deltat to get a weight.
    """
    value = used_misfit_window.deltat
    weight = cal_cos_weight_ops(value, value1, value2)
    return weight

##################################################################################
# weights used considering all the windowed traces
##################################################################################


def cal_geographical_weight(stations_mapper, used_net_sta_list):
    """
    cal_geographical_weight: calculate geographycal weighting for the given stations.
        + stations_mapper: stations_mapper[net_sta]=(lon,lat)
        + used_net_sta_list: a list of net_sta used.
    """
    # * build up the lon,lat matrix in the order of used_net_sta_list
    matrix_size = len(used_net_sta_list)
    lat1 = np.zeros((matrix_size, matrix_size))
    lon1 = np.zeros((matrix_size, matrix_size))
    lat2 = np.zeros((matrix_size, matrix_size))
    lon2 = np.zeros((matrix_size, matrix_size))
    for irow in range(matrix_size):
        for icolumn in range(matrix_size):
            lat1[irow, icolumn] = stations_mapper[used_net_sta_list[irow]][1]
            lon1[irow, icolumn] = stations_mapper[used_net_sta_list[irow]][0]
            lat2[irow, icolumn] = stations_mapper[used_net_sta_list[icolumn]][1]
            lon2[irow, icolumn] = stations_mapper[used_net_sta_list[icolumn]][0]
    # calculate the distance matrix
    distance_matrix = locations2degrees(lat1, lon1, lat2, lon2)

    # * find Delta0 and get geographical weighting for each net_sta
    dref = np.arange(0.5, 8.5, 0.5)
    ratio_list = []
    for iref in range(len(dref)):
        wt = np.zeros(matrix_size)
        for i in range(matrix_size):
            wt[i] = 1.0/np.sum(np.exp(-(distance_matrix[i]/dref[iref])**2))
        ratio = np.max(wt)/np.min(wt)
        ratio_list.append(ratio)
    # we should find the nearest 1/3 max ratio
    ratio_list = np.array(ratio_list)
    max_ratio = np.max(ratio_list)
    candiate_ratio = max_ratio/3
    candiate_id = np.argmin(
        np.abs(ratio_list[:len(ratio_list)//2]-candiate_ratio))
    used_ref = dref[candiate_id]
    used_wt_dict = {}
    for i in range(matrix_size):
        net_sta = used_net_sta_list[i]
        used_wt_dict[net_sta] = 1.0 / \
            np.sum(np.exp(-(distance_matrix[i]/used_ref)**2))
    return used_wt_dict


def cal_category_weight(number_in_category):
    if (number_in_category == 0):
        return 0
    else:
        return 1.0/number_in_category
