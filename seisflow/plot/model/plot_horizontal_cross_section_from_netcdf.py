"""
plot_horizontal_cross_section_from_netcdf.py: plot the horizontal cross section from the netcdf model file.
"""
import cartopy
import cartopy.crs as ccrs
import click
import matplotlib.pyplot as plt
import numpy as np
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from scipy.interpolate import RegularGridInterpolator
from scipy.io import netcdf


def extract_data(f, depth, parameter):
    depth_pos = np.where(f.variables["depth"][:] == depth)
    if (len(depth_pos) != 1):
        raise Exception("no such depth in the netcdf file.")
    depth_pos = depth_pos[0][0]
    data = f.variables[parameter][:, :, depth_pos].copy()
    # we usually use a large value to represent nan.
    data[data > 100] = np.nan
    mesh_lon, mesh_lat = np.meshgrid(
        f.variables["longitude"][:], f.variables["latitude"][:], indexing="ij")
    return mesh_lon, mesh_lat, data


def plot_h(mesh_lon, mesh_lat, data, parameter, depth, vmin, vmax, region):
    plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    plt.pcolormesh(mesh_lon, mesh_lat, data,
                   transform=ccrs.PlateCarree(), vmin=vmin, vmax=vmax, cmap=plt.cm.jet_r)  # pylint: disable=no-member
    ax.coastlines()
    # input format lon1/lat1/lon2/lat2
    lon1, lat1, lon2, lat2 = map(float, region.split("/"))
    ax.set_extent([lon1, lon2, lat1, lat2], crs=ccrs.PlateCarree())
    ax.set_xticks(np.arange(lon1, lon2, 10), crs=ccrs.PlateCarree())
    ax.set_yticks(np.arange(lat1, lat2, 10), crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    # in general disable it, since the map boundary excludes some part of China.
    ax.add_feature(cartopy.feature.BORDERS)
    plt.colorbar(orientation='horizontal', fraction=0.046, pad=0.04)
    plt.grid()
    plt.title(f"{parameter} at {depth}km")

    plt.show()


@click.command()
@click.option('--netcdf_file', required=True, type=str, help="the netcdf file")
@click.option('--parameter', required=True, type=str, help="the parameter to plot")
@click.option('--depth', required=True, type=float, help="the depth to extract data (km)")
@click.option('--vmin', required=True, type=float, help="the min limit for colorbar")
@click.option('--vmax', required=True, type=float, help="the max limit for colorbar")
@click.option('--region', required=True, type=str, help="the region to plot, lon1/lat1/lon2/lat2")
def main(netcdf_file, parameter, depth, vmin, vmax, region):
    f = netcdf.netcdf_file(netcdf_file, 'r')
    mesh_lon, mesh_lat, data = extract_data(f, depth, parameter)
    plot_h(mesh_lon, mesh_lat, data, parameter, depth, vmin, vmax, region)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
