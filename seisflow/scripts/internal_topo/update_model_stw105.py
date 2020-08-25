"""
update_model_stw105.py: update the d660 topo and stw105 ppm based on the kernel text file.
"""
import click
import numpy as np
from scipy.interpolate import griddata

from .generate_gaussian_stw105 import (generate_mapper, generate_ppm,
                                       read_stw105, write_to_netcdf)


@click.command()
@click.option('--lat_range', required=True, type=str, help="latmin,latmax should be larger than the simulation region")
@click.option('--lon_range', required=True, type=str, help="lonmin,lonmax should be larger than the simulation region")
@click.option('--kernel_path', required=True, type=str, help="the generated kernel text file from julia")
@click.option('--output_path_d660', required=True, type=str, help="the external mesh file output path")
@click.option('--output_path_ppm', required=True, type=str, help="the external ppm file output path")
def main(lat_range, lon_range, kernel_path, output_path_d660, output_path_ppm):
    lat_min, lat_max = map(float, lat_range.split(","))
    lon1 = np.arange(-180, 180.2, 0.2)
    lat1 = np.arange(lat_min, lat_max + 0.2, 0.2)
    lon2, lat2 = np.meshgrid(lon1, lat1, indexing="ij")

    # * read in the model and kernel
    raw_data = np.loadtxt(kernel_path)
    # * interp
    input_kernel = griddata(raw_data[:, :2], raw_data[:, 3],
                            (lon2, lat2), method='linear')
    input_depth = griddata(
        raw_data[:, :2], raw_data[:, 2], (lon2, lat2), method='linear')
    # for the nan part of the inputs, we use the predefined values
    input_kernel[np.isnan(input_kernel)] = 0
    input_depth[np.isnan(input_depth)] = 650

    # * we should normalize the kernel so the maximum absolute value can be corresponding to 5km
    max_abs_pos = np.unravel_index(
        np.argmax(np.abs(input_kernel)), input_kernel.shape)
    ratio = 5 / np.abs(input_kernel[max_abs_pos])
    # update the kernel based on the ratio
    output_kernel = ratio * input_kernel

    # * add the kernel to the depth
    output_depth = input_depth + output_kernel

    # * now we write to the topo file, need to change lons from 0 to 360
    towrite_depth = np.zeros_like(output_depth)
    towrite_depth[:900, :] = output_depth[901:, :]
    towrite_depth[900:, :] = output_depth[:901, :]
    towrite_depth = towrite_depth-650
    towrite_lon1 = np.arange(0, 360.2, 0.2)
    towrite_lat1 = lat1
    with open(output_path_d660, "w") as f:
        f.write(f"{towrite_depth.shape[0]} {towrite_depth.shape[1]} \n")
        for ilat in range(len(towrite_lat1)):
            for ilon in range(len(towrite_lon1)):
                f.write(
                    f"{towrite_lon1[ilon]:.2f} {towrite_lat1[ilat]:.2f} 0.0000 {towrite_depth[ilon,ilat]:.4f} \n")

    # * now we are ready to generate the ppm file based in the d660 file
    latmin, latmax = map(float, lat_range.split(","))
    lonmin, lonmax = map(float, lon_range.split(","))
    # based on lonmin and lonmax, we extract the data from output_depth
    lonmin_index = int((lonmin + 180) * 5)
    lonmax_index = int((lonmax + 180) * 5)
    d650 = output_depth[lonmin_index:(lonmax_index+1), :]

    # * now we generate the ppm file
    stw105_iso = read_stw105()
    mapper_all, mapper_depression, mapper_uplifting = generate_mapper(
        d650, stw105_iso)
    lon1, lat1, dep1, vs3d, vp3d, rho3d = generate_ppm(
        d650, (latmin, latmax), (lonmin, lonmax), mapper_all, mapper_depression, mapper_uplifting)
    write_to_netcdf(lon1, lat1, dep1, vs3d, vp3d, rho3d, output_path_ppm)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
