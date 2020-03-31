"""
generate the source mask list.
"""
from glob import glob
from os.path import join

import numpy as np


@click.command()
@click.option('--base_directory', required=True, type=str, help="the base simulation directory")
@click.option('--output_path', required=True, type=str, help="the output path for the source mask list")
def main(base_directory, output_path):
    all_event_dirs = glob(join(base_directory, "*"))
    source_mask_list = []
    for each_event_dir in all_event_dirs:
        source_mask_vtk_path = join(
            each_event_dir, "OUTPUT_FILES", "source.vtk")
        source_mask = np.loadtxt(source_mask_vtk_path, skiprows=5)
        source_mask_list.append(source_mask)
    source_mask_list = np.array(source_mask_list)
    np.savetxt(output_path, source_mask_list)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
