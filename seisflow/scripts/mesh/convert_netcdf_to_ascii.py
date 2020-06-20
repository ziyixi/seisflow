"""
convert_netcdf_to_ascii.py: convert the netcdf files to the ascii format.
"""
from os.path import join

import click
import numpy as np
from mpi4py import MPI
from netCDF4 import Dataset

comm = MPI.COMM_WORLD  # pylint: disable=c-extension-no-member
size = comm.Get_size()
rank = comm.Get_rank()


@click.command()
@click.option('--netcdf_file', required=True, type=str)
@click.option('--output_dir', required=True, type=str)
@click.option('--models', required=False, default="vs,vp,vsv,vsh,vpv,vph,eta,rho", type=str)
def main(netcdf_file, output_dir, models):
    models = models.split(",")
    models_this_rank = np.array_split(models, size)[rank]
    with Dataset(netcdf_file, 'r') as f:
        lat = f.variables["latitude"][:]
        lon = f.variables["longitude"][:]
        dep = f.variables["depth"][:]
        for each_parameter in models_this_rank:
            data = f.variables[each_parameter][:]
            nlon, nlat, ndep = data.shape
            with open(join(output_dir, each_parameter), "w") as g:
                for ilon in range(nlon):
                    for ilat in range(nlat):
                        for idep in range(ndep):
                            g.write(
                                f"{lon[ilon]:.2f} {lat[ilat]:.2f} {dep[idep]:.1f} {data[ilon,ilat,idep]:.5f}\n")


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
