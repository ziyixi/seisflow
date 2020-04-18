"""
mpi_simplify_data.py: simplify data asdf after processing.
"""
from glob import glob
from os.path import join

import numpy as np
from mpi4py import MPI

from ...tasks.asdf.simplify_data_asdf import simplify_data_asdf

comm = MPI.COMM_WORLD  # pylint: disable=c-extension-no-member
size = comm.Get_size()
rank = comm.Get_rank()


def get_files_this_rank(asdf_directory):
    all_files = sorted(glob(join(asdf_directory, "*h5")))
    all_files_this_rank = np.array_split(all_files, size)[rank]
    return all_files_this_rank


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--asdf_directory', required=True, type=str, help="asdf files directory.")
    def main(asdf_directory):
        all_files_this_rank = get_files_this_rank(asdf_directory)
        for each_file in all_files_this_rank:
            simplify_data_asdf(each_file)
    main()  # pylint: disable=no-value-for-parameter
