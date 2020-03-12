"""
plot_vertical_cross_section_from_netcdf.py: plot the vertical cross section from the netcdf file.
"""
from scipy.io import netcdf
from scipy.interpolate import RegularGridInterpolator
import numpy as np
import matplotlib.pyplot as plt


def generate_vcs_mesh(lon1, lat1, lon2, lat2, dep1, dep2, rh, rdep, theta_label):
    if(theta_label == "lat"):
        if(lat1 < lat2):
            lat_list = np.arange(lat1, lat2+rh, rh)
        else:
            lat_list = np.arange(lat1, lat2-rh, -rh)
        ntheta = len(lat_list)
        theta_list = lat_list
        lon_list = np.linspace(lon1, lon2, ntheta)
    elif(theta_label == "lon"):
        if(lon1 < lon2):
            lon_list = np.arange(lon1, lon2+rh, rh)
        else:
            lon_list = np.arange(lon1, lon2-rh, -rh)
        ntheta = len(lon_list)
        theta_list = lon_list
        lat_list = np.linspace(lat1, lat2, ntheta)
    theta_range = range(ntheta)
    dep_list = np.arange(dep1, dep2+rdep, rdep)
    ndep = len(dep_list)
    dep_range = range(ndep)
    array_to_interpolate = np.zeros((ntheta, ndep, 3))
    for itheta in theta_range:
        for idep in dep_range:
            array_to_interpolate[itheta, idep, :] = [
                lon_list[itheta], lat_list[itheta], dep_list[idep]]
    # generate mesh
    mesh_theta, mesh_dep = np.meshgrid(
        theta_list, dep_list, indexing="ij")
    return mesh_theta, mesh_dep, array_to_interpolate


def get_interp_function(netcdf_data, parameter, method="nearest"):
    interpolating_function = RegularGridInterpolator(
        (netcdf_data.variables["longitude"][:], netcdf_data.variables["latitude"][:], netcdf_data.variables["depth"][:]), netcdf_data.variables[parameter][:], method=method)
    return interpolating_function


def extract_data_v(netcdf_data, lon1, lat1, lon2, lat2, dep1, dep2, rh, rdep, theta_label, parameter, method="nearest"):
    mesh_theta, mesh_dep, array_to_interpolate = generate_vcs_mesh(
        lon1, lat1, lon2, lat2, dep1, dep2, rh, rdep, theta_label)
    interpolating_function = get_interp_function(
        netcdf_data, parameter, method=method)
    plot_values = interpolating_function(array_to_interpolate)
    return mesh_theta, mesh_dep, plot_values


def plot_v(lat1, lat2, lon1, lon2, theta_label, mesh_theta, mesh_dep, plot_values):
    fig = plt.figure(figsize=(20, 14))
    ax = fig.add_subplot(111, polar=True)
    if(theta_label == "lat"):
        lat2, lat1 = min(lat1, lat2), max(lat1, lat2)
        ax.set_thetamin(lat2)
        ax.set_thetamax(lat1)
        ax.set_theta_zero_location("N", offset=-(lat1+lat2)/2)
    else:
        lon2, lon1 = min(lon1, lon2), max(lon1, lon2)
        ax.set_thetamin(lon2)
        ax.set_thetamax(lon1)
        ax.set_theta_zero_location("N", offset=-(lon1+lon2)/2)
    ax.set_rorigin(6371)
    ax.set_rlim(bottom=800, top=0)
    # ax.scatter(40,300,s=100)

    contourf_ = ax.pcolormesh(np.deg2rad(
        mesh_theta), mesh_dep, plot_values, cmap=plt.cm.seismic_r, vmin=-0.06, vmax=0.06)
    # plt.colorbar(contourf_)
    plt.colorbar(contourf_, orientation='horizontal', fraction=0.046, pad=-0.4)
    ax.set_ylabel('Depth(km)', size=20)

    plt.show()
