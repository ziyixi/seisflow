"""
mpi_generate_perturbed_sync_given_step_length.py: generate the perturbed sync given step length, raw sync and perturbed sync (setting).
"""
from glob import glob
from os.path import basename, join

import click
import numpy as np
from mpi4py import MPI

from ...tasks.line_search.line_search_structure import get_perturbed_vir_sync
from ...utils.asdf_io import VirAsdf

comm = MPI.COMM_WORLD  # pylint: disable=c-extension-no-member
size = comm.Get_size()
rank = comm.Get_rank()


def get_raw_per_asdf_this_rank(all_raw_asdf, all_per_asdf):
    assert len(all_raw_asdf) == len(all_per_asdf)
    all_raw_asdf_this_rank = np.array_split(all_raw_asdf, size)[rank]
    all_per_asdf_this_rank = np.array_split(all_per_asdf, size)[rank]
    return all_raw_asdf_this_rank, all_per_asdf_this_rank


@click.command()
@click.option('--raw_directory', required=True, type=str, help="raw asdf directory")
@click.option('--perturbed_directory', required=True, type=str, help="perturbed asdf directory (based on step length in setting)")
@click.option('--output_directory', required=True, type=str, help="output directory")
@click.option('--step_length', required=True, type=float, help="asdf of step length that is supposed to generate")
def main(raw_directory, perturbed_directory, output_directory, step_length):
    all_raw_asdf = sorted(glob(join(raw_directory, "*h5")))
    all_per_asdf = sorted(glob(join(perturbed_directory, "*h5")))
    all_raw_asdf_this_rank, all_per_asdf_this_rank = get_raw_per_asdf_this_rank(
        all_raw_asdf, all_per_asdf)
    assert len(all_raw_asdf_this_rank) == len(all_per_asdf_this_rank)
    for each_raw, each_per in zip(all_raw_asdf_this_rank, all_per_asdf_this_rank):
        assert basename(each_raw) == basename(each_per)
        thebase = basename(each_raw)
        output_path = join(output_directory, thebase)
        virasdf_raw = VirAsdf()
        virasdf_raw.read_asdf(each_raw)
        virasdf_perturbed = VirAsdf()
        virasdf_perturbed.read_asdf(each_per)
        virasdf_vir_perturbed = get_perturbed_vir_sync(
            virasdf_raw, virasdf_perturbed, step_length)
        virasdf_vir_perturbed.write_asdf(output_path)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
