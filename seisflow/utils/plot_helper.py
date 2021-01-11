""" 
Functions to help plotting figures.
"""
from functools import partial

import numpy as np
import pyproj
from obspy.geodetics import degrees2kilometers
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize


def gmt_lon_project(lon_start, lon_end, lon0, lat0, az):
    """
    Get the points along the great circle line with the even spacing of the longitude.
    """
    g = pyproj.Geod(ellps='WGS84')
    # * firstly we find the starting points and the ending points by finding the root.

    def find_lat_raw(lat, lon, lat0, lon0):
        # get the azimuth between lon0,lat0 and lon,lat
        used_lat = lat
        used_lon = np.ones_like(used_lat)*lon
        used_lat0 = np.ones_like(used_lat)*lat0
        used_lon0 = np.ones_like(used_lat)*lon0
        azs = np.array(g.inv(used_lon0, used_lat0, used_lon, used_lat)[0])
        # * azs might not be in the range of 0~360, here we convert it
        azs = azs % 360
        # * replace nan with the nearest non-nan
        diff = azs - az
        # * **2 to find 0
        return diff**2

    # * generate the result lon list
    result_lons = np.linspace(lon_start, lon_end, 1000)

    # * for each lon, find the lat that meets the az requirement
    result_lats = []
    for each_lon in result_lons:
        find_lat = partial(find_lat_raw, lon=each_lon, lat0=lat0, lon0=lon0)
        sol = minimize(find_lat, 0, bounds=[(-90, 90)])
        result_lats.append(sol.x[0])
    result_lats = np.array(result_lats)
    # * fix the issue when lon0==lon_start
    if(lon0 == lon_start):
        result_lats[0] = lat0

    return result_lons, result_lats


def gmt_lat_project(lat_start, lat_end, lon0, lat0, az):
    """
    Get the points along the great circle line with the even spacing of the latitude.
    """
    g = pyproj.Geod(ellps='WGS84')
    # * firstly we find the starting points and the ending points by finding the root.

    def find_lon_raw(lon, lat, lat0, lon0):
        # get the azimuth between lon0,lat0 and lon,lat
        used_lon = lon
        used_lat = np.ones_like(used_lon)*lat
        used_lat0 = np.ones_like(used_lon)*lat0
        used_lon0 = np.ones_like(used_lon)*lon0
        azs = np.array(g.inv(used_lon0, used_lat0, used_lon, used_lat)[0])
        # * azs might not be in the range of 0~360, here we convert it
        azs = azs % 360
        # * replace nan with the nearest non-nan
        diff = azs - az
        # * **2 to find 0
        return diff**2

    # * generate the result lon list
    result_lats = np.linspace(lat_start, lat_end, 1000)

    # * for each lon, find the lat that meets the az requirement
    result_lons = []
    for each_lat in result_lats:
        find_lon = partial(find_lon_raw, lat=each_lat, lat0=lat0, lon0=lon0)
        sol = minimize(find_lon, 0, bounds=[(-360, 360)])
        result_lons.append(sol.x[0])
    result_lons = np.array(result_lons)

    return result_lons, result_lats


def gmt_dist_project(dist_start, dist_end, lon0, lat0, az):
    """
    Get the points along the great circle line with the even spacing of the point distances.
    """
    g = pyproj.Geod(ellps='WGS84')
    # * note when we are talking distance in degree of Earth, it should refer to the spherical Earth or meaningless
    # firstly we convert the distance to meters
    dist_start_meter = degrees2kilometers(dist_start)*1000
    dist_end_meter = degrees2kilometers(dist_end)*1000

    # get the starting and endding points
    lon_start, lat_start, _ = g.fwd(lon0, lat0, az, dist_start_meter)
    lon_end, lat_end, _ = g.fwd(lon0, lat0, az, dist_end_meter)

    result = g.npts(lon_start, lat_start, lon_end, lat_end, 1000)
    result = np.array(result)

    return result[:, 0], result[:, 1]


def model_interp(to_interp_data, lons, lats, deps):
    """
    Give an xarray model, interp it based on the given lats, lons, deps and construct a new xarray dataset.
    mainly used to generate the vertical cross-sections
    """
    # * len(lons) should be the same as len(lats)
    profile_list = []
    for idep in range(len(deps)):
        for ilon in range(len(lons)):
            profile_list.append([lons[ilon], lats[ilon], deps[idep]])
    model_interpolating_function = RegularGridInterpolator(
        (to_interp_data.longitude.data, to_interp_data.latitude.data, to_interp_data.depth.data), to_interp_data.data)
    interp_result = model_interpolating_function(profile_list)
    cross_section = np.zeros((len(lons), len(deps)))

    icount = 0
    for idep in range(len(deps)):
        for ilon in range(len(lons)):
            cross_section[ilon, idep] = interp_result[icount]
            icount += 1

    # cross_section_xarray = xr.DataArray(cross_section, dims=(
    #     'h', "v"), coords={'h': lons, "v": deps})

    return cross_section


def topo_interp(to_interp_data, lons, lats):
    """
    Give the xarray topography model, interp the elevation line along the given (lons,lats) pair. 
    """
    profile_list = []
    for ilon in range(len(lons)):
        profile_list.append([lons[ilon], lats[ilon]])
    # the names and the transverse might be adjusted, this is the gmt format
    grd_interpolating_function = RegularGridInterpolator(
        (to_interp_data.lon.data, to_interp_data.lat.data), to_interp_data.data.T)

    grd_interp_result = grd_interpolating_function(profile_list)

    # * return the 1d array
    return grd_interp_result
