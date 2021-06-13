"""
add_isotropic_velocities.py: add vp and vs velocity for the absolute velocity netcdf files.
"""
import click
import numpy as np
import xarray as xr


@click.command()
@click.option('--target_netcdf', required=True, type=str)
def main(target_netcdf):
    with xr.open_dataset(target_netcdf) as data:
        data.load()
    vsv = data["vsv"].data
    vsh = data["vsh"].data
    vpv = data["vpv"].data
    vph = data["vph"].data
    # reference  doi: 10.1111/j.1365-246X.2006.03100.x
    vs = np.sqrt((2 * vsv ** 2 + vsh ** 2) / 3)
    vp = np.sqrt((vpv ** 2 + 4 * vph ** 2) / 5)
    data["vs"] = (('longitude', 'latitude', 'depth'), vs)
    data["vp"] = (('longitude', 'latitude', 'depth'), vp)
    data.to_netcdf(target_netcdf)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
