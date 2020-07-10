"""
add_isotropic_velocities.py: add vp and vs velocity for the absolute velocity netcdf files.
"""
import click
from netCDF4 import Dataset


@click.command()
@click.option('--target_netcdf', required=True, type=str)
@click.option('--output_netcdf', required=True, type=str)
@click.option('--models', required=False, default="vs,vp,vpv,vph,vsv,vsh,eta,rho", type=str)
def main(target_netcdf, output_netcdf, models):
    models = models.split(",")
    with Dataset(target_netcdf, 'r') as target:
        with Dataset(output_netcdf, 'w') as f:
            f.createDimension(
                'depth', len(target.variables["depth"][::2]))
            f.createDimension(
                'latitude', target.variables["latitude"][:].shape[0])
            f.createDimension(
                'longitude', target.variables["longitude"][:].shape[0])
            longitude = f.createVariable('longitude', 'f8', ('longitude',))
            longitude[:] = target.variables["longitude"][:]
            latitude = f.createVariable('latitude', 'f8', ('latitude',))
            latitude[:] = target.variables["latitude"][:]
            depth = f.createVariable('depth', 'f8', ('depth',))
            depth[:] = target.variables["depth"][::2]
            for index_parameter in range(len(models)):
                parameter_name = models[index_parameter]
                target_array = target.variables[parameter_name][:].copy()
                parameter_array = target_array[:, :, ::2]
                netcdf_var = f.createVariable(
                    parameter_name, 'f8', ('longitude', 'latitude', 'depth'), zlib=True)
                netcdf_var[:] = parameter_array


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
