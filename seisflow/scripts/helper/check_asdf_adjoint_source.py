"""
check_asdf_adjoint_source.py: check if the adjoint sources have abnormal values in a directory
"""
from glob import glob
from os.path import basename, join

import click
import numpy as np
import pyasdf
from mpi4py import MPI

comm = MPI.COMM_WORLD  # pylint: disable=c-extension-no-member
size = comm.Get_size()
rank = comm.Get_rank()


def get_asdf_files_this_rank(adjoint_directory):
    """
    get asdf files to read this rank.
    """
    all_files = sorted(glob(join(adjoint_directory, "*h5")))
    adjoint_files_this_rank = np.array_split(all_files, size)[rank]
    return adjoint_files_this_rank


def kernel_check_file(adjoint_file_path):
    with pyasdf.ASDFDataSet(adjoint_file_path, mode="r", mpi=False) as ds:
        all_adjoint_ids = ds.auxiliary_data.AdjointSources.list()
        for each_adjoint_id in all_adjoint_ids:
            adjoint_data = ds.auxiliary_data.AdjointSources[each_adjoint_id].data[:]
            status = np.isfinite(adjoint_data).all()
            if (not status):
                print(basename(adjoint_file_path), each_adjoint_id)


@click.command()
@click.option('--adjoint_directory', required=True, type=str, help="the adjoint files directory")
def main(adjoint_directory):
    adjoint_files_this_rank = get_asdf_files_this_rank(adjoint_directory)
    for each_adjoint_file in adjoint_files_this_rank:
        kernel_check_file(each_adjoint_file)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
