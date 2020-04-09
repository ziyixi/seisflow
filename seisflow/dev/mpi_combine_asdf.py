"""
mpi_combine_asdf.py: combine two asdf files together.
"""
from os.path import join, basename, isfile
from glob import glob
from ..asdf.combine_asdf_files import combine_asdf
import numpy as np
from mpi4py import MPI
import click

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


def get_asdf_pairs(base_directory, append_directory, output_directory):
    """
    The file names in the twp directory should be the same.
    """
    all_base_paths = sorted(glob(join(base_directory, "*")))
    all_append_paths = sorted(glob(join(append_directory, "*")))
    assert len(all_base_paths) == len(all_append_paths)
    pair_list = []
    for each_base_path, each_append_path in zip(all_base_paths, all_append_paths):
        assert basename(each_base_path) == basename(each_append_path)
        output_path = join(output_directory, basename(each_base_path))
        pair_list.append((each_base_path, each_append_path, output_path))
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
        combine_asdf(each_pair[0], each_pair[1], each_pair[2])


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
