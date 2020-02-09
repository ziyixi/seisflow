"""
ln_new_model_to_gll.py: remove the new flag and ln the generated gll files.
"""
from glob import glob
from os.path import basename, join

import click
import sh


@click.command()
@click.option('--new_flag_dir', required=True, type=str, help="the directory containning the models with the flag new")
@click.option('--output_dir', required=True, type=str, help="the output directory")
def main(new_flag_dir, output_dir):
    """
    main: the main program to ln the kernels.
    """
    all_kernels = sorted(glob(join(new_flag_dir, "*bin")))
    sh.mkdir("-p", output_dir)
    for each_smoothed_kernel_path in all_kernels:
        base_fname = basename(each_smoothed_kernel_path)
        base_fname_splitted = base_fname.split("_")
        if base_fname_splitted[-1] != "new.bin":
            continue
        new_fname = "_".join(base_fname_splitted[:-1]) + ".bin"
        to_ln_path = join(output_dir, new_fname)
        sh.ln("-s", each_smoothed_kernel_path, to_ln_path)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
