"""
convert_txt_output_2_netcdf.py: convert the output ppm to the netcdf file.
"""
import numpy as np
import tqdm
import click
import numba
from os.path import join
from numba import float64
from scipy.io import netcdf


@numba.njit
def pos_mapper(lon, lat, dep, slon, slat, sdep, dlon, dlat, ddep):
    i = int(round((lon-slon)/dlon, 0))
    j = int(round((lat-slat)/dlat, 0))
    k = int(round((dep-sdep)/ddep, 0))
    return [i, j, k]


@click.command()
@click.option('--base_dir', required=True, type=str)
@click.option('--region', required=True, type=str, help="minlon/maxlon/minlat/maxlat/mindep/maxdep")
@click.option('--npts', required=True, type=str, help="nlon/nlat/ndep")
@click.option('--nproc', required=True, type=int)
@click.option('--parameters', required=True, type=str)
@click.option('--out_path', required=True, type=str)
@click.option('--history', required=True, type=str, help="information to write in the file")
def main(base_dir, region, npts, nproc, parameters, out_path, history):
    minlon, maxlon, minlat, maxlat, mindep, maxdep = [
        float(item) for item in region.split("/")]
    nlon, nlat, ndep = [int(item) for item in npts.split("/")]
    dlon = (maxlon-minlon)/(nlon-1)
    dlat = (maxlat-minlat)/(nlat-1)
    ddep = (maxdep-mindep)/(ndep-1)
    parameters_list = parameters.split(",")
    save_arrays_list = [np.zeros((nlon, nlat, ndep))
                        for i in range(len(parameters_list))]

    # for each file, we load the file, and map it to the parameter numpy array
    for iproc in tqdm.tqdm(range(nproc)):
        fname = join(base_dir, str(iproc))
        data_iproc = np.loadtxt(fname)
        # for each parameter array, map the value
        for index_point in range(data_iproc.shape[0]):
            for index_parameter in range(len(parameters_list)):
                i, j, k = pos_mapper(
                    data_iproc[index_point, 0], data_iproc[index_point, 1], data_iproc[index_point, 2], minlon, minlat, mindep, dlon, dlat, ddep)
                save_arrays_list[index_parameter][i,
                                                  j, k] = data_iproc[index_point, index_parameter+3]

    # save to netcdf file
    with netcdf.netcdf_file(out_path, 'w') as f:
        f.history = history
        f.createDimension('depth', ndep)
        f.createDimension('latitude', nlat)
        f.createDimension('longitude', nlon)
        longitude = f.createVariable('longitude', 'f8', ('longitude',))
        longitude[:] = np.linspace(minlon, maxlon, nlon)
        latitude = f.createVariable('latitude', 'f8', ('latitude',))
        latitude[:] = np.linspace(minlat, maxlat, nlat)
        depth = f.createVariable('depth', 'f8', ('depth',))
        depth[:] = np.linspace(mindep, maxdep, ndep)
        for index_parameter in range(len(parameters_list)):
            parameter_array = save_arrays_list[index_parameter]
            parameter_name = parameters_list[index_parameter]
            netcdf_var = f.createVariable(
                parameter_name, 'f8', ('longitude', 'latitude', 'depth'))
            netcdf_var[:] = parameter_array


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
