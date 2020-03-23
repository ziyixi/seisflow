"""
plot_vertical_cross_section_from_netcdf.py: plot the vertical cross section from the netcdf file.
"""
import click
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.io import netcdf


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


def plot_v(lat1, lat2, lon1, lon2, dep2, theta_label, mesh_theta, mesh_dep, plot_values, vmin, vmax):
    fig = plt.figure()
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
    ax.set_rorigin(6000)
    ax.set_rlim(bottom=dep2, top=0)
    plt.subplots_adjust(bottom=0.2)
    print(np.min(plot_values), np.max(plot_values))
    contourf_ = ax.pcolormesh(np.deg2rad(
        mesh_theta), mesh_dep, (plot_values), cmap=plt.cm.jet_r, vmin=vmin, vmax=vmax)  # pylint: disable=no-member
    plt.colorbar(contourf_, orientation='horizontal', fraction=0.046, pad=0.0)

    plt.show()


@click.command()
@click.option('--netcdf_file', required=True, type=str, help="the netcdf file")
@click.option('--parameter', required=True, type=str, help="the parameter to plot")
@click.option('--vmin', required=True, type=float, help="the min colorbar threshold")
@click.option('--vmax', required=True, type=float, help="the max colorbar threshold")
@click.option('--region', required=True, type=str, help="plot region, lon1/lat1/lon2/lat2/dep1/dep2")
@click.option('--rh', required=True, type=float, help="the horizontal resolution")
@click.option('--rdep', required=True, type=float, help="the vertical resolution")
@click.option('--theta_label', required=True, type=str, help="can be lat or lon")
def main(netcdf_file, parameter, vmin, vmax, region, rh, rdep, theta_label):
    with netcdf.netcdf_file(netcdf_file, 'r') as f:
        interpolating_function = get_interp_function(
            f, parameter, method="linear")
        lon1, lat1, lon2, lat2, dep1, dep2 = map(float, region.split("/"))
        mesh_theta, mesh_dep, array_to_interpolate = generate_vcs_mesh(
            lon1, lat1, lon2, lat2, dep1, dep2, rh, rdep, theta_label)
        plot_values = interpolating_function(array_to_interpolate)
        plot_v(lat1, lat2, lon1, lon2, dep2, theta_label,
               mesh_theta, mesh_dep, plot_values, vmin, vmax)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
