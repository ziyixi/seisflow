"""
make_perturbed_netcdf.py: make the perturbed netcdf based on two nc files.
"""
from scipy.io import netcdf
import click


@click.command()
@click.option('--target_netcdf', required=True, type=str)
@click.option('--reference_netcdf', required=True, type=str)
@click.option('--output_netcdf', required=True, type=str)
@click.option('--models', required=False, default="vpv,vph,vsv,vsh,eta,rho", type=str)
def main(target_netcdf, reference_netcdf, output_netcdf, models):
    models = models.split(",")
    with netcdf.netcdf_file(target_netcdf, 'r') as target:
        with netcdf.netcdf_file(reference_netcdf, 'r') as reference:
            with netcdf.netcdf_file(output_netcdf, 'w') as f:
                history = f"perturbation for {target_netcdf} respect to {reference_netcdf} "
                f.history = history
                f.createDimension(
                    'depth', target.variables["depth"][:].shape[0])
                f.createDimension(
                    'latitude', target.variables["latitude"][:].shape[0])
                f.createDimension(
                    'longitude', target.variables["longitude"][:].shape[0])
                longitude = f.createVariable('longitude', 'f8', ('longitude',))
                longitude[:] = target.variables["longitude"][:]
                latitude = f.createVariable('latitude', 'f8', ('latitude',))
                latitude[:] = target.variables["latitude"][:]
                depth = f.createVariable('depth', 'f8', ('depth',))
                depth[:] = target.variables["depth"][:]
                for index_parameter in range(len(models)):
                    parameter_name = models[index_parameter]
                    # get parameter_array
                    parameter_array = (
                        target.variables[parameter_name][:] - reference.variables[parameter_name][:]) / reference.variables[parameter_name][:]
                    parameter_array[target.variables[parameter_name]
                                    [:] > 9e6] = 9999999.
                    netcdf_var = f.createVariable(
                        parameter_name, 'f8', ('longitude', 'latitude', 'depth'))
                    netcdf_var[:] = parameter_array


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
