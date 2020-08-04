"""
plot_horizontal_cross_section_from_netcdf.py: plot the horizontal cross section from the netcdf model file.
"""
import cartopy
import cartopy.crs as ccrs
import click
import matplotlib.pyplot as plt
import numpy as np
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from netCDF4 import Dataset


def extract_data(f, depth, parameter):
    depth_pos = np.where(f.variables["depth"][:] == depth)
    if (len(depth_pos) != 1):
        raise Exception("no such depth in the netcdf file.")
    depth_pos = depth_pos[0][0]
    data_all = f.variables[parameter][:, :, :].copy()
    # data_all=data_all.transpose([2,1,0])
    data = data_all[:, :, depth_pos].copy()
    # we usually use a large value to represent nan.
    data[data > 9e6] = np.nan
    data_all[data_all > 9e6] = np.nan
    print(np.nanmin(data_all), np.nanmax(data_all))
    mesh_lon, mesh_lat = np.meshgrid(
        f.variables["longitude"][:], f.variables["latitude"][:], indexing="ij")
    return mesh_lon, mesh_lat, data


def plot_h(mesh_lon, mesh_lat, data, parameter, depth, vmin, vmax, region, scale, percentage, abs):  # pylint: disable=unused-argument
    plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    print(np.nanmin(data), np.nanmax(data))
    if (scale):
        min_absdata = np.abs(np.nanmin(data))
        max_absdata = np.abs(np.nanmax(data))
        range_val = np.max([min_absdata, max_absdata])
        vmin = -range_val
        vmax = range_val
    if (abs):
        data = np.abs(data)
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
    colorbar = plt.colorbar(orientation='horizontal', fraction=0.046,
                            pad=0.04, extend="neither")
    colorbar.set_ticks([round(i, 2) for i in colorbar.get_ticks().tolist()])
    colorbar.set_ticklabels(
        [str(i) for i in colorbar.get_ticks().tolist()])
    # colorbar.set_label(label=f"${{dlnV_{{{parameter[1:]}}}}}$", size=15)
    plt.grid()
    plt.text(0.1, 0.9, f"{parameter}\n{depth}km", ha='center',
             va='center', transform=ax.transAxes, fontsize=20)
    ax.xaxis.set_tick_params(labelsize=15)
    ax.yaxis.set_tick_params(labelsize=15)
    # plt.title(f"{parameter} at {depth}km")

    plt.show()


@click.command()
@click.option('--netcdf_file', required=True, type=str, help="the netcdf file")
@click.option('--parameter', required=True, type=str, help="the parameter to plot")
@click.option('--depth', required=True, type=float, help="the depth to extract data (km)")
@click.option('--vmin', required=False, default=0, type=float, help="the min limit for colorbar")
@click.option('--vmax', required=False, default=0, type=float, help="the max limit for colorbar")
@click.option('--region', required=True, type=str, help="the region to plot, lon1/lat1/lon2/lat2")
@click.option('--scale/--no-scale', default=False, required=False, help="if scale the range based on the maximum value")
@click.option('--percentage/--no-percentage', default=False, required=False, help="if use percentage in colorbar")
@click.option('--abs/--no-abs', default=False, required=False, help="if plot the absolute value")
def main(netcdf_file, parameter, depth, vmin, vmax, region, scale, percentage, abs):
    f = Dataset(netcdf_file, 'r')
    mesh_lon, mesh_lat, data = extract_data(f, depth, parameter)
    if (percentage):
        data = data*100
    plot_h(mesh_lon, mesh_lat, data, parameter,
           depth, vmin, vmax, region, scale, percentage, abs)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
