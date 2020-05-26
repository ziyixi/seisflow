"""
convert_netcdf_to_ascii.py: convert the netcdf files to the ascii format.
"""
from scipy.io import netcdf
import click
from os.path import join


@click.command()
@click.option('--netcdf_file', required=True, type=str)
@click.option('--output_dir', required=True, type=str)
@click.option('--models', required=False, default="vs,vp,vsv,vsh,vpv,vph,eta,rho", type=str)
def main(netcdf_file, output_dir, models):
    models = models.split(",")
    with netcdf.netcdf_file(netcdf_file, 'r') as f:
        lat = f.variables["latitude"][:]
        lon = f.variables["longitude"][:]
        dep = f.variables["depth"][:]
        for each_parameter in models:
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
