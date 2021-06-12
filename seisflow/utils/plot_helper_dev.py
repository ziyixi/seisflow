""" 
Functions to help plotting figures.
"""
import numpy as np
import pyproj
from scipy.spatial import KDTree
import tempfile
from obspy.geodetics import locations2degrees


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

def gmt_lon_as_dist(startlon, startlat, endlon, endlat,a_list=None,g_interval=None,npts=1001):
    g = pyproj.Geod(ellps='WGS84')
    gcarc=locations2degrees(startlat,startlon,endlat,endlon)
    if startlon>endlon:
        startlon,endlon=endlon,startlon
        startlat,endlat=endlat,startlat
    test_points = np.array(g.npts(startlon, startlat, endlon, endlat, (npts-1)*10+1))
    tree = KDTree(test_points[:, 0].reshape(test_points.shape[0], -1))
    # * generate the gmt custome axis file when use lon to plot the cross-section
    num_list=[]
    type_list=[]
    annote_list=[]
    # before first a
    for each_lon in np.arange(a_list[0]-g_interval,startlon,-1*g_interval)[::-1]:
        num_list.append(each_lon)
        type_list.append("f")
        annote_list.append("")
    # start from first a
    for ia in range(len(a_list)-1):
        starta=a_list[ia]
        enda=a_list[ia+1]
        # first a
        num_list.append(starta)
        type_list.append("a")
        annote_list.append(f"{starta}")
        for each_g in np.arange(starta+g_interval,enda,g_interval):
            num_list.append(each_g)
            type_list.append("f")
            annote_list.append("")
    # last a
    num_list.append(a_list[-1])
    type_list.append("a")
    annote_list.append(f"{a_list[-1]}")
    for each_g in np.arange(a_list[-1]+g_interval,endlon,g_interval):
        num_list.append(each_g)
        type_list.append("f")
        annote_list.append("")

    # convert num_list to actual dist_list
    num_list=np.array(num_list)
    _, pos = tree.query(num_list.reshape(len(num_list), -1))
    dist_list=pos/((npts-1)*10)*gcarc

    # write lists to a temp file
    tmp = tempfile.NamedTemporaryFile()
    with open(tmp.name, 'w') as f:
        for index in range(len(num_list)):
            f.write(f"{dist_list[index]}  {type_list[index]}  {annote_list[index]} \n")
    return tmp