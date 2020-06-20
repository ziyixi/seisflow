"""
add_isotropic_velocities.py: add vp and vs velocity for the absolute velocity netcdf files.
"""
import click
import numpy as np
from netCDF4 import Dataset


@click.command()
@click.option('--target_netcdf', required=True, type=str)
def main(target_netcdf):
    with Dataset(target_netcdf, 'a') as f:
        # * vsv,vsh,vpv,vph must exist in the file, and we generate vs and vp.
        vsv = f.variables["vsv"][:]
        vsh = f.variables["vsh"][:]
        vpv = f.variables["vpv"][:]
        vph = f.variables["vph"][:]
        # calculate
        # reference  doi: 10.1111/j.1365-246X.2006.03100.x
        vs = np.sqrt((2 * vsv ** 2 + vsh ** 2) / 3)
        vp = np.sqrt((vpv ** 2 + 4 * vph ** 2) / 5)
        # write vp and vs to the netcdf file
        netcdf_vp = f.createVariable(
            "vp", 'f8', ('longitude', 'latitude', 'depth'))
        netcdf_vs = f.createVariable(
            "vs", 'f8', ('longitude', 'latitude', 'depth'))
        netcdf_vp[:] = vp
        netcdf_vs[:] = vs


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
