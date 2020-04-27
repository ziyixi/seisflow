"""
convert_txt_output_2_netcdf.py: convert the output ppm to the netcdf file.
"""
from os.path import join

import click
import numba
import numpy as np
import tqdm
from mpi4py import MPI
from scipy.io import netcdf

comm = MPI.COMM_WORLD  # pylint: disable=c-extension-no-member
size = comm.Get_size()
rank = comm.Get_rank()


@numba.njit
def pos_mapper(lon, lat, dep, slon, slat, sdep, dlon, dlat, ddep):
    i = int(round((lon-slon)/dlon, 0))
    j = int(round((lat-slat)/dlat, 0))
    k = int(round((dep-sdep)/ddep, 0))
    return [i, j, k]


def get_pos_mapper_collection(base_dir, iproc, minlon, minlat, mindep, dlon, dlat, ddep, parameters_list):
    """
    For each proc, make the mapper list to use later. (i,j,k,index_point,index_parameter)
    """
    fname = join(base_dir, str(iproc))
    data_iproc = np.loadtxt(fname)
    result_list = []
    # for each parameter array, map the value
    for index_point in range(data_iproc.shape[0]):
        for index_parameter in range(len(parameters_list)):
            i, j, k = pos_mapper(
                data_iproc[index_point, 0], data_iproc[index_point, 1], data_iproc[index_point, 2], minlon, minlat, mindep, dlon, dlat, ddep)
            value = data_iproc[index_point, index_parameter+3]
            result_list.append([i, j, k, value, index_parameter])
    result_list = np.array(result_list)
    return result_list


def get_iproc_this_rank(all_iproc):
    return np.array_split(all_iproc, size)[rank]


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

    # save_arrays_list = [np.zeros((nlon, nlat, ndep))
    #                     for i in range(len(parameters_list))]
    # for each file, we load the file, and map it to the parameter numpy array
    # for iproc in tqdm.tqdm(range(nproc)):
    #     fname = join(base_dir, str(iproc))
    #     data_iproc = np.loadtxt(fname)
    #     # for each parameter array, map the value
    #     for index_point in range(data_iproc.shape[0]):
    #         for index_parameter in range(len(parameters_list)):
    #             i, j, k = pos_mapper(
    #                 data_iproc[index_point, 0], data_iproc[index_point, 1], data_iproc[index_point, 2], minlon, minlat, mindep, dlon, dlat, ddep)
    #             save_arrays_list[index_parameter][i,
    #                                             j, k] = data_iproc[index_point, index_parameter+3]

    # * we get the iproc for this rank
    iproc_this_rank = get_iproc_this_rank(range(nproc))
    pos_mapper_collection_this_rank = []
    if(rank == 0):
        pbar = tqdm.tqdm(total=len(iproc_this_rank))
    for iproc in iproc_this_rank:
        pos_mapper_iproc_rank = get_pos_mapper_collection(
            base_dir, iproc, minlon, minlat, mindep, dlon, dlat, ddep, parameters_list)
        pos_mapper_collection_this_rank.append(pos_mapper_iproc_rank)
        if(rank == 0):
            pbar.update(1)
    pos_mapper_collection_this_rank = np.vstack(
        pos_mapper_collection_this_rank)
    # * get the maximum array row length
    all_row_lengths = comm.gather(
        pos_mapper_collection_this_rank.shape[0], root=0)
    max_row_length = np.max(all_row_lengths)
    if(rank == 0):
        pos_mapper_collection = np.zeros((size, max_row_length, 5))
        pos_mapper_collection[:] = np.nan
    else:
        pos_mapper_collection = None
    # * collect all the pos_mapper_collection
    comm.Gather(
        pos_mapper_collection_this_rank, pos_mapper_collection, root=0)
    if(rank == 0):
        pos_mapper_collection = np.vstack(pos_mapper_collection)
        # * now we map to save_arrays_list
        save_arrays_list = [np.zeros((nlon, nlat, ndep))
                            for i in range(len(parameters_list))]
        for row in pos_mapper_collection:
            i, j, k, value, index_parameter = row
            if (np.isnan(i)):
                continue
            save_arrays_list[index_parameter][int(i), int(j), int(k)] = value

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
