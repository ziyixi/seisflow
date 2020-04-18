"""
mpi_process_sync_structure_inversion.py: process the sync in series, but parallel for multiple events.
"""
from glob import glob
from os.path import basename, join

import click
import numpy as np
from mpi4py import MPI

from ...tasks.process.process_sync_single_st import process_sync_single_trace
from ...utils.asdf_io import VirAsdf

comm = MPI.COMM_WORLD  # pylint: disable=c-extension-no-member
size = comm.Get_size()
rank = comm.Get_rank()


def get_asdf_this_rank(sync_directory):
    all_files = sorted(glob(join(sync_directory, "*h5")))
    asdf_this_rank = np.array_split(all_files, size)[rank]
    return asdf_this_rank


def process_single_asdf(asdf_file_name, output_directory, waveform_length, taper_tmin_tmax, min_period, max_period, sampling_rate):
    # * load asdf file
    asdf = VirAsdf()
    asdf.read_asdf(asdf_file_name, usecopy=True)
    # * loop for each trace, process the sync
    waveforms_list = asdf.get_waveforms_list()
    event = asdf.get_events()[0]
    origin = event.preferred_origin() or event.origins[0]
    event_time = origin.time
    waveforms = asdf.get_waveforms()
    tag = "preprocessed_%is_to_%is" % (
        int(min_period), int(max_period))
    for net_sta in waveforms_list:
        st = waveforms[net_sta]["st"]
        process_sync_single_trace(st, event_time, waveform_length,
                                  taper_tmin_tmax, min_period, max_period, sampling_rate)
    gcmtid = basename(asdf_file_name).split(".")[0]
    output_fname = f"{gcmtid}.{tag}.h5"
    asdf.write_asdf(join(output_directory, output_fname), tag)


@click.command()
@click.option('--sync_directory', required=True, type=str, help="the raw sync directory")
@click.option('--output_directory', required=True, type=str, help="the processed sync directory")
@click.option('--periods', required=True, type=str, help="min periods in filtering: minp1,maxp1/minp2,maxp2/...")
@click.option('--waveform_length', required=True, type=str, help="the length of the waveform to cut")
@click.option('--sampling_rate', required=True, type=int, help="the sampling rate to use")
@click.option('--taper_tmin_tmaxs', required=True, type=str, help="the taper time bands: minp1,maxp1/minp2,maxp2/...")
def main(sync_directory, output_directory, periods, waveform_length, sampling_rate, taper_tmin_tmaxs):
    periods_split = periods.split("/")
    taper_tmin_tmaxs_split = taper_tmin_tmaxs.split("/")
    assert len(periods_split) == len(taper_tmin_tmaxs_split)
    asdf_this_rank = get_asdf_this_rank(sync_directory)
    for each_period, each_taper in zip(periods_split, taper_tmin_tmaxs_split):
        for each_asdf in asdf_this_rank:
            min_period, max_period = map(float, each_period.split(","))
            process_single_asdf(each_asdf, output_directory, waveform_length,
                                each_taper, min_period, max_period, sampling_rate)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
