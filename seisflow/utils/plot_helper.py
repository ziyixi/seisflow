""" 
Functions to help plotting figures.
"""
from functools import partial

import numpy as np
import pyproj
from obspy.geodetics import degrees2kilometers
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
