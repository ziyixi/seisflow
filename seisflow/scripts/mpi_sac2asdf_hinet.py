"""
mpi_sac2asdf_hinet.py: convert the downloaded files from hinet downloader to the asdf format.
"""
from glob import glob
from os.path import basename, isdir, join

import numpy as np
import sh
from mpi4py import MPI

from .sac2asdf_hinet import Sac2Asdf

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


def unzip_single_event(event_directory):
    """
    unzip the sac files in a directory
    """
    current_path = str(sh.pwd())[:-1]  # pylint: disable=not-callable
    sh.cd(event_directory)
    sh.tar("xvf", "PZ.tar.gz")
    sh.cd(current_path)


def rm_single_event(event_directory):
    """
    rm the extracted directory
    """
    current_path = str(sh.pwd())[:-1]  # pylint: disable=not-callable
    sh.cd(event_directory)
    sh.rm("-rf", "PZ")
    sh.cd(current_path)


def get_directories_this_rank(base_directory):
    """
    get directories to convert in this rank.
    """
    all_directories = sorted(glob(join(base_directory, "*")))
    all_directories = [item for item in all_directories if isdir(item)]
    directories_this_rank = np.array_split(all_directories, size)[rank]
    return directories_this_rank


def convert_this_rank(directories_this_rank, cmt_directory, output_directory):
    """
    convert the sac files to the asdf format in this rank.
    """
    for event_directory in directories_this_rank:
        unzip_single_event(event_directory)
        gcmtid = basename(event_directory)
        cmt_path = join(cmt_directory, gcmtid)
        output_path = join(output_directory, gcmtid)
        run_script = Sac2Asdf(
            event_directory, cmt_path, output_path)
        run_script.run()
        rm_single_event(event_directory)


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--base_directory', required=True, type=str, help="the directory containing all the event data directories")
    @click.option('--cmt_directory', required=True, type=str, help="the directory for cmt files")
    @click.option('--output_directory', required=True, type=str, help="the output asdf directory")
    def main(base_directory, cmt_directory, output_directory):
        """
        convert all the sac files for the Hi-net to the asdf format.
        """
        directories_this_rank = get_directories_this_rank(base_directory)
        convert_this_rank(directories_this_rank,
                          cmt_directory, output_directory)

    main()  # pylint: disable=no-value-for-parameter
