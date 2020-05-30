"""
mpi_combine_asdf.py: combine two asdf files together.
"""
from glob import glob
from os.path import basename, join, isfile

import click
import numpy as np
from mpi4py import MPI

from ...tasks.asdf.combine_asdf_files import combine_asdf
import sh

comm = MPI.COMM_WORLD  # pylint: disable=c-extension-no-member
size = comm.Get_size()
rank = comm.Get_rank()


def get_asdf_pairs(base_directory, append_directory, output_directory):
    """
    The file names in the twp directory should be the same.
    """
    # * we will always base on the number of the base directory
    all_base_paths = sorted(glob(join(base_directory, "*")))
    pair_list = []
    for each_base_path in all_base_paths:
        thebase = basename(each_base_path)
        each_append_path = join(append_directory, thebase)
        output_path = join(output_directory, thebase)
        if (isfile(each_append_path)):
            pair_list.append((each_base_path, each_append_path, output_path))
        else:
            pair_list.append((each_base_path, output_path))
    return pair_list


def get_asdf_pairs_this_rank(pair_list):
    pair_list_this_rank = np.array_split(pair_list, size)[rank]
    return pair_list_this_rank


@click.command()
@click.option('--base_directory', required=True, type=str, help="the base asdf directory")
@click.option('--append_directory', required=True, type=str, help="the append asdf directory")
@click.option('--output_directory', required=True, type=str, help="the output asdf directory")
def main(base_directory, append_directory, output_directory):
    pair_list = get_asdf_pairs(
        base_directory, append_directory, output_directory)
    pair_list_this_rank = get_asdf_pairs_this_rank(pair_list)
    for each_pair in pair_list_this_rank:
        if (len(each_pair) == 2):
            sh.cp(each_pair[0], each_pair[1])
        elif(len(each_pair) == 3):
            combine_asdf(each_pair[0], each_pair[1], each_pair[2])
        else:
            raise Exception("len(each_pair) is not 2 or 3")


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
