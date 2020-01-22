"""
mpi_generate_waveform_comparision_source_inversion.py: compare the waveform during the source inversion massively.
"""
from glob import glob
from os.path import basename, join

import click
from mpi4py import MPI

from ...utils.get_path import get_asdf_fnames
from .compare_two_asdf_with_windows import kernel
import numpy as np

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


@click.command()
@click.option('--cmts_directory', required=True, type=str)
@click.option('--obs_directory', required=True, type=str)
@click.option('--syn_directory', required=True, type=str)
@click.option('--output_directory', required=True, type=str)
@click.option('--info_dir', required=True, type=str)
@click.option('--misfit_windows_directory', required=True, type=str)
@click.option('--min_periods', required=True, type=str)
@click.option('--max_periods', required=True, type=str)
@click.option('--azimuth_width', required=True, type=int)
@click.option('--waves_perpage', required=True, type=int)
@click.option('--trace_length', required=True, type=int)
@click.option('--snr', required=True, type=float)
@click.option('--cc', required=True, type=float)
@click.option('--deltat', required=True, type=float)
def main(cmts_directory, obs_directory, syn_directory, output_directory, info_dir, misfit_windows_directory,
         min_periods, max_periods, azimuth_width, waves_perpage, trace_length, snr, cc, deltat):
    # get gcmtids
    cmts_paths = glob(join(cmts_directory, "*"))
    all_gcmtids = [basename(item) for item in cmts_paths]
    all_gcmtids_this_rank = np.array_split(all_gcmtids, size)[rank]
    for gcmtid in all_gcmtids_this_rank:
        data_asdf_body_path, sync_asdf_body_path, data_asdf_surface_path, sync_asdf_surface_path = get_asdf_fnames(
            gcmtid, min_periods, max_periods, obs_directory, syn_directory)
        # * firstly we generate pdf file for the body wave band
        body_min_period, surface_min_period = min_periods.split(",")
        body_max_period, surface_max_period = max_periods.split(",")
        output_pdf_path = join(
            output_directory, f"{gcmtid}.preprocessed_{float(body_min_period):.0f}s_to_{float(body_max_period):.0f}s.pdf")
        misfit_windows_path = join(misfit_windows_directory, f"{gcmtid}.pkl")
        kernel(data_asdf_body_path, sync_asdf_body_path,
               azimuth_width, output_pdf_path, waves_perpage, trace_length, info_dir, misfit_windows_path, snr, cc, deltat, False)


if __name__ == "__main__":
    main()
