"""
ln_smoothed_kernel_to_input_dir: remove the flag smooth and link the smoothed kernels to the directory INPUT_GRADIENT.
"""
from glob import glob
from os.path import basename, join

import click
import sh


@click.command()
@click.option('--smooth_dir', required=True, type=str, help="the directory containning the smoothed kernels")
@click.option('--kernel_process_directory', required=True, type=str, help="the specfem directory to process the kernels")
def main(smooth_dir, kernel_process_directory):
    """
    main: the main program to ln the kernels.
    """
    all_kernels = sorted(glob(join(smooth_dir, "*bin")))
    sh.mkdir("-p", join(kernel_process_directory, "INPUT_GRADIENT"))
    for each_smoothed_kernel_path in all_kernels:
        base_fname = basename(each_smoothed_kernel_path)
        base_fname_splitted = base_fname.split("_")
        # if base_fname_splitted[-1] != "smooth.bin":
        #     raise Exception("check if the bin file is the smoothed kernels!")
        # new_fname = "_".join(base_fname_splitted[:-1]) + ".bin"
        # * note since the following script need smooth flag, we just do a soft link here.
        new_fname = "_".join(base_fname_splitted)
        to_ln_path = join(kernel_process_directory,
                          "INPUT_GRADIENT", new_fname)
        try:
            sh.unlink(to_ln_path)
        except sh.ErrorReturnCode_1:
            pass
        sh.ln("-s", each_smoothed_kernel_path, to_ln_path)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
