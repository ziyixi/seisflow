"""
make_perturbed_cmtsolution.py: make perturbed cmt solution using the frechet information.
"""
from glob import glob
from os.path import basename, join

import numpy as np
import obspy

from ...tasks.cmtsolution.make_perturbed_cmtsolution import add_src_frechet
from ...utils.save_files import save_cmtsolution
from ...setting import MAX_DXS_RATIO

# * test is passed for this script on 01/07/2020


def load_src_frechet(src_frechet_path):
    return np.loadtxt(src_frechet_path)


def load_cmtsolution(cmtsolution_path):
    cmtsolution = obspy.read_events(cmtsolution_path)[0]
    return cmtsolution


def kernel(src_frechet_directory, cmtsolution_directory, output_directory, max_dxs_ratio, fix_location=False, fix_focal=False):
    # we assume the gcmtids in src_frechet_directory and cmtsolution_directory are the same
    all_files_frechet = sorted(glob(join(src_frechet_directory, "*")))
    all_gcmtids_frechet = [basename(item).split(
        ".")[0] for item in all_files_frechet]
    all_files_cmtsolution = sorted(glob(join(cmtsolution_directory, "*")))
    all_gcmtids_cmtsolution = [basename(item).split(
        ".")[0] for item in all_files_cmtsolution]
    all_gcmtids = sorted(set(all_gcmtids_frechet) &
                         set(all_gcmtids_cmtsolution))

    for each_gcmtid in all_gcmtids:
        # get paths
        src_frechet_path = join(src_frechet_directory, each_gcmtid)
        cmtsolution_path = join(cmtsolution_directory, each_gcmtid)
        output_path = join(output_directory, each_gcmtid)
        # run
        src_frechet = load_src_frechet(src_frechet_path)
        cmtsolution = load_cmtsolution(cmtsolution_path)
        cmtsolution_new = add_src_frechet(
            src_frechet, cmtsolution, max_dxs_ratio, fix_location=fix_location, fix_focal=fix_focal)
        save_cmtsolution(output_path, cmtsolution_new)


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--src_frechet_directory', required=True, type=str, help="src_frechet files directory")
    @click.option('--cmtsolution_directory', required=True, type=str, help="cmtsolution files directory")
    @click.option('--output_directory', required=True, type=str, help="output CMTSOLUTION directory")
    @click.option('--fix_location/--no-fix_location', default=False, required=False, help="if fix location for the perturbed cmtsolution")
    @click.option('--fix_focal/--no-fix_focal', default=False, required=False, help="if fix focal mechanism for the perturbed cmtsolution")
    def main(src_frechet_directory, cmtsolution_directory, output_directory, fix_location, fix_focal):
        max_dxs_ratio = MAX_DXS_RATIO
        kernel(src_frechet_directory, cmtsolution_directory,
               output_directory, max_dxs_ratio, fix_location=fix_location, fix_focal=fix_focal)
    main()  # pylint: disable=no-value-for-parameter
