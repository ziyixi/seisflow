"""
make_perturbed_netcdf.py: make the perturbed netcdf based on two nc files.
"""
import xarray as xr
import click
import numpy as np


@click.command()
@click.option('--target_netcdf', required=True, type=str)
@click.option('--reference_netcdf', required=True, type=str)
@click.option('--output_netcdf', required=True, type=str)
@click.option('--models', required=False, default="vpv,vph,vsv,vsh,eta,rho", type=str)
def main(target_netcdf, reference_netcdf, output_netcdf, models):
    models = models.split(",")
    with xr.open_dataset(target_netcdf) as target:
        with xr.open_dataset(reference_netcdf) as reference:
            for each_parameter in models:
                target[each_parameter].data = target[each_parameter].data / \
                    reference[each_parameter].data-1
            target.to_netcdf(output_netcdf)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
