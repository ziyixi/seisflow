"""
generate_gaussian_stw105.py: generate the ppm model for stw105 with perturbed d660.
"""
import click
import numpy as np
from netCDF4 import Dataset  # pylint: disable=no-name-in-module
from scipy import interpolate

from ... import get_data
from .generate_gaussian_d660 import generate_650_topo


def read_stw105():
    """
    read the stw105 raw file and convert to the isotropic case.
    """
    stw105 = np.loadtxt(get_data("stw105.txt"), dtype=np.float)
    stw105_iso = np.zeros((stw105.shape[0], 4))
    for index, row in enumerate(stw105[::-1, :]):
        depth = (6371000-row[0])/1000
        rho = row[1]/1000
        vpv, vsv = row[2]/1000, row[3]/1000
        vph, vsh = row[6]/1000, row[7]/1000
        vs = np.sqrt((2 * vsv ** 2 + vsh ** 2) / 3)
        vp = np.sqrt((vpv ** 2 + 4 * vph ** 2) / 5)
        stw105_iso[index, :] = [depth, vs, vp, rho]
    return stw105_iso


def generate_mapper(d650, stw105_iso):
    """
    generate the map from the depth to the model, considering extrapolation from above or from below of d660.
    """
    # depths
    d650_max = np.int(np.trunc(np.max(d650))+20)
    d650_min = np.int(np.trunc(np.min(d650)) - 20)

    # * generate mapper for all.
    vs_all = interpolate.interp1d(stw105_iso[:, 0], stw105_iso[:, 1])
    vp_all = interpolate.interp1d(stw105_iso[:, 0], stw105_iso[:, 2])
    rho_all = interpolate.interp1d(stw105_iso[:, 0], stw105_iso[:, 3])
    mapper_all = {}
    # we consider the problem down to the outer core
    for thedepth in range(0, 2893):
        mapper_all[thedepth] = {
            "vs": vs_all(thedepth),
            "vp": vp_all(thedepth),
            "rho": rho_all(thedepth)
        }

    # * generate mapper of depression, from d650_min to 650km
    stw105_iso_depression = stw105_iso[np.where(
        ((650 + d650_min) < stw105_iso[:, 0]) & (stw105_iso[:, 0] < 650))]
    vs_depression = interpolate.interp1d(
        stw105_iso_depression[:, 0], stw105_iso_depression[:, 1], fill_value="extrapolate")
    vp_depression = interpolate.interp1d(
        stw105_iso_depression[:, 0], stw105_iso_depression[:, 2], fill_value="extrapolate")
    rho_depression = interpolate.interp1d(
        stw105_iso_depression[:, 0], stw105_iso_depression[:, 3], fill_value="extrapolate")
    mapper_depression = {}
    for thedepth in range(650, 650 + d650_max+1):
        mapper_depression[thedepth] = {
            "vs": vs_depression(thedepth),
            "vp": vp_depression(thedepth),
            "rho": rho_depression(thedepth)
        }

    # * generate mapper of uplifting, from 650km to d650_max
    stw105_iso_uplifting = stw105_iso[np.where(
        ((650 + d650_max) > stw105_iso[:, 0]) & (stw105_iso[:, 0] > 650))]
    vs_uplifting = interpolate.interp1d(
        stw105_iso_uplifting[:, 0], stw105_iso_uplifting[:, 1], fill_value="extrapolate")
    vp_uplifting = interpolate.interp1d(
        stw105_iso_uplifting[:, 0], stw105_iso_uplifting[:, 2], fill_value="extrapolate")
    rho_uplifting = interpolate.interp1d(
        stw105_iso_uplifting[:, 0], stw105_iso_uplifting[:, 3], fill_value="extrapolate")
    mapper_uplifting = {}
    for thedepth in range(650 + d650_min, 651):
        mapper_uplifting[thedepth] = {
            "vs": vs_uplifting(thedepth),
            "vp": vp_uplifting(thedepth),
            "rho": rho_uplifting(thedepth)
        }

    return mapper_all, mapper_depression, mapper_uplifting


def generate_ppm(d650, lats, lons, mapper_all, mapper_depression, mapper_uplifting):
    """
    generate the ppm 3d matrix, lons, lats and deps.
    """
    lon1 = np.arange(lons[0], lons[1] + 0.2, 0.2)
    lat1 = np.arange(lats[0], lats[1] + 0.2, 0.2)
    dep1 = np.arange(0, 2893)
    vs3d = np.zeros((2893, len(lat1), len(lon1)))
    vp3d = np.zeros((2893, len(lat1), len(lon1)))
    rho3d = np.zeros((2893, len(lat1), len(lon1)))

    d650_max = np.int(np.trunc(np.max(d650))+20)
    d650_min = np.int(np.trunc(np.min(d650)) - 20)

    # * for each depth, we consider the position of the point, and use appropriate mapper to get the model value
    for idep in range(2893):
        thedepth = idep
        if (thedepth < 650 + d650_min):
            vs3d[idep, :, :] = mapper_all[thedepth]["vs"]
            vp3d[idep, :, :] = mapper_all[thedepth]["vp"]
            rho3d[idep, :, :] = mapper_all[thedepth]["rho"]
        elif (thedepth > 650 + d650_max):
            vs3d[idep, :, :] = mapper_all[thedepth]["vs"]
            vp3d[idep, :, :] = mapper_all[thedepth]["vp"]
            rho3d[idep, :, :] = mapper_all[thedepth]["rho"]
        else:
            # * we have to consider the problem point by point
            for ilat in range(len(lat1)):
                for ilon in range(len(lon1)):
                    d650_depth = 650 + d650[ilon, ilat]
                    if (d650_depth <= thedepth < 650):
                        # * uplifting
                        vs3d[idep, ilat, ilon] = mapper_uplifting[thedepth]["vs"]
                        vp3d[idep, ilat, ilon] = mapper_uplifting[thedepth]["vp"]
                        rho3d[idep, ilat, ilon] = mapper_uplifting[thedepth]["rho"]
                    elif(650 < thedepth <= d650_depth):
                        # * depression, we have the value at 650km
                        vs3d[idep, ilat, ilon] = mapper_depression[thedepth]["vs"]
                        vp3d[idep, ilat, ilon] = mapper_depression[thedepth]["vp"]
                        rho3d[idep, ilat, ilon] = mapper_depression[thedepth]["rho"]
                    elif (thedepth == 650):
                        if (d650_depth < 650):
                            vs3d[idep, ilat, ilon] = mapper_uplifting[thedepth]["vs"]
                            vp3d[idep, ilat, ilon] = mapper_uplifting[thedepth]["vp"]
                            rho3d[idep, ilat,
                                  ilon] = mapper_uplifting[thedepth]["rho"]
                        else:
                            vs3d[idep, ilat, ilon] = mapper_depression[thedepth]["vs"]
                            vp3d[idep, ilat, ilon] = mapper_depression[thedepth]["vp"]
                            rho3d[idep, ilat,
                                  ilon] = mapper_depression[thedepth]["rho"]
                    else:
                        vs3d[idep, ilat, ilon] = mapper_all[thedepth]["vs"]
                        vp3d[idep, ilat, ilon] = mapper_all[thedepth]["vp"]
                        rho3d[idep, ilat, ilon] = mapper_all[thedepth]["rho"]
    return lon1, lat1, dep1, vs3d, vp3d, rho3d


def write_to_netcdf(lon1, lat1, dep1, vs3d, vp3d, rho3d, output_path):
    """
    write the model files to a netcdf file
    """
    with Dataset(output_path, 'w') as f:
        f.createDimension('depth', len(dep1)+1)
        f.createDimension('latitude', len(lat1))
        f.createDimension('longitude', len(lon1))
        longitude = f.createVariable('longitude', 'f4', ('longitude',))
        longitude[:] = lon1
        latitude = f.createVariable('latitude', 'f4', ('latitude',))
        latitude[:] = lat1
        depth = f.createVariable('depth', 'f4', ('depth',))
        # set value above 0 to avoid issue in interpolating
        depth[:] = np.hstack([dep1[0]-dep1[1], dep1])
        vs = f.createVariable(
            "vs", 'f4', ('depth', 'latitude', 'longitude'), zlib=True, fill_value=999999)
        vs3d_input = np.zeros([vs3d.shape[0]+1, vs3d.shape[1], vs3d.shape[2]])
        vs3d_input[1:, :, :] = vs3d
        vs3d_input[0, :, :] = vs3d[0, :, :]
        vs[:] = vs3d_input
        vp = f.createVariable(
            "vp", 'f4', ('depth', 'latitude', 'longitude'), zlib=True, fill_value=999999)
        vp3d_input = np.zeros([vp3d.shape[0]+1, vp3d.shape[1], vp3d.shape[2]])
        vp3d_input[1:, :, :] = vp3d
        vp3d_input[0, :, :] = vp3d[0, :, :]
        vp[:] = vp3d_input
        rho = f.createVariable(
            "rho", 'f4', ('depth', 'latitude', 'longitude'), zlib=True, fill_value=999999)
        rho3d_input = np.zeros(
            [rho3d.shape[0]+1, rho3d.shape[1], rho3d.shape[2]])
        rho3d_input[1:, :, :] = rho3d
        rho3d_input[0, :, :] = rho3d[0, :, :]
        rho[:] = rho3d_input


@click.command()
@click.option('--lats', required=True, type=str, help="latmin,latmax should be larger than the simulation region")
@click.option('--lons', required=True, type=str, help="lonmin,lonmax should be larger than the simulation region")
@click.option('--dep0', required=True, type=float, help="the depression(+) or uplifting(-) depth")
@click.option('--lon0', required=True, type=float, help="the longitude for the depression/uplifting")
@click.option('--lat0', required=True, type=float, help="the latitude for the depression/uplifting")
@click.option('--radius', required=True, type=float, help="the radius for the depression/uplifting")
@click.option('--output_path', required=True, type=str, help="the output path for the depression/uplifting ppm file")
def main(lons, lats, dep0, lon0, lat0, radius, output_path):
    # * firstly we generate the d650
    latmin, latmax = map(float, lats.split(","))
    lonmin, lonmax = map(float, lons.split(","))
    d650 = generate_650_topo(latmin, latmax, dep0, lon0,
                             lat0, radius, None, lonmin, lonmax)
    # * stw105
    stw105_iso = read_stw105()
    # * get the model and write to a netcdf file
    mapper_all, mapper_depression, mapper_uplifting = generate_mapper(
        d650, stw105_iso)
    lon1, lat1, dep1, vs3d, vp3d, rho3d = generate_ppm(
        d650, (latmin, latmax), (lonmin, lonmax), mapper_all, mapper_depression, mapper_uplifting)
    write_to_netcdf(lon1, lat1, dep1, vs3d, vp3d, rho3d, output_path)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
