"""
generate the source mask list.
"""
from glob import glob
from os.path import join

import click
import numpy as np


@click.command()
@click.option('--base_directory', required=True, type=str, help="the base simulation directory")
@click.option('--source_mask_km', required=True, type=float, help="radius for source mask in km")
@click.option('--receiver_mask_km', required=True, type=float, help="radius for receiver mask in km")
def main(base_directory, source_mask_km, receiver_mask_km):
    all_event_dirs = glob(join(base_directory, "*"))
    for each_event_dir in all_event_dirs:
        source_mask_vtk_path = join(
            each_event_dir, "OUTPUT_FILES", "source.vtk")
        receiver_mask_vtk_path = join(
            each_event_dir, "OUTPUT_FILES", "receiver.vtk")
        source_mask = np.loadtxt(source_mask_vtk_path, skiprows=5)
        source_mask = np.array([source_mask])
        source_mask_append = np.ones(
            (np.shape(source_mask)[0], 1)) * source_mask_km
        source_mask = np.hstack([source_mask, source_mask_append])
        receiver_mask = np.loadtxt(receiver_mask_vtk_path, skiprows=5)
        receiver_mask_append = np.ones(
            (receiver_mask.shape[0], 1)) * receiver_mask_km
        receiver_mask = np.hstack([receiver_mask, receiver_mask_append])
        mask_list = np.vstack([source_mask, receiver_mask])
        output_path = join(
            each_event_dir, "OUTPUT_FILES", "mask.xyz")
        np.savetxt(output_path, mask_list, fmt="%.6E %.6E %.6E %.6E")


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
