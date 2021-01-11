""" 
Functions to help plotting figures.
"""
import numpy as np
import pyproj
from scipy.spatial import KDTree


def gmt_project(startlon, startlat, endlon, endlat, thetype, npts=1001):
    g = pyproj.Geod(ellps='WGS84')
    if(thetype == "dist"):
        result = g.npts(startlon, startlat, endlon, endlat, npts)
        result = np.array(result)
        # result_lons, result_lats
        return result[:, 0], result[:, 1]
    elif(thetype == "lon"):
        # we divide 10 more times points using npts and find nearest lon
        test_points = g.npts(startlon, startlat, endlon, endlat, (npts-1)*10+1)
        test_points = np.array(test_points)
        tree = KDTree(test_points[:, 0].reshape(test_points.shape[0], -1))
        # evenly distributed lons
        result_lons = np.linspace(startlon, endlon, npts)
        _, pos = tree.query(result_lons.reshape(npts, -1))
        result_lats = test_points[:, 1][pos]
        return result_lons, result_lats
    elif(thetype == "lat"):
        # we divide 10 more times points using npts and find nearest lat
        test_points = g.npts(startlon, startlat, endlon, endlat, (npts-1)*10+1)
        test_points = np.array(test_points)
        tree = KDTree(test_points[:, 1].reshape(test_points.shape[0], -1))
        # evenly distributed lons
        result_lats = np.linspace(startlat, endlat, npts)
        _, pos = tree.query(result_lats.reshape(npts, -1))
        result_lons = test_points[:, 0][pos]
        return result_lons, result_lats
    else:
        raise Exception("not supported type")
