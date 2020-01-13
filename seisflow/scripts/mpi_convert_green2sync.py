"""
mpi_convert_green2sync.py: convert the green function asdf to the normal convolved asdf in a directory.
"""
from ..tasks.process.convert_green2processed import conv_process_single_virasdf
from ..utils.asdf_io import VirAsdf
from mpi4py import MPI
from glob import glob
from os.path import join, basename
import numpy as np
import obspy

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


def get_paths_this_rank(cmts_directory, green_directory, conv_directory):
    all_files = glob(join(green_directory, "*h5"))
    asdf_paths_this_rank = np.array_split(all_files, size)[rank]
    cmt_paths_this_rank = []
    conv_paths_this_rank = []
    for each_asdf_path_this_rank in asdf_paths_this_rank:
        gcmtid = basename(each_asdf_path_this_rank).split(".")[0]
        cmt_paths_this_rank.append(join(cmts_directory, gcmtid))
        conv_paths_this_rank.append(join(conv_directory, gcmtid))
    return asdf_paths_this_rank, cmt_paths_this_rank, conv_paths_this_rank


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--cmts_directory', required=True, type=str, help="the directory of the cmtsolution")
    @click.option('--green_directory', required=True, type=str, help="the directory of the green function sync")
    @click.option('--conv_directory', required=True, type=str, help="the directory to output the convolution result")
    @click.option('--waveform_length', required=True, type=str, help="the length of the waveform to cut")
    @click.option('--taper_tmin_tmaxs', required=True, type=str, help="the taper time bands: minp1,maxp1/minp2,maxp2/...")
    @click.option('--periods', required=True, type=str, help="min periods in filtering: minp1,maxp1/minp2,maxp2/...")
    @click.option('--sampling_rate', required=True, type=int, help="the sampling rate to use")
    def main(cmts_directory, green_directory, conv_directory,
             waveform_length, taper_tmin_tmaxs, periods, sampling_rate):
        # firstly we handle the input parameters
        taper_tmin_tmax_list = taper_tmin_tmaxs.split("/")
        periods_list = periods.split("/")
        asdf_paths_this_rank, cmt_paths_this_rank, conv_paths_this_rank = get_paths_this_rank(
            cmts_directory, green_directory, conv_directory)
        for each_taper_tmin_tmax, each_period in zip(taper_tmin_tmax_list, periods_list):
            taper_tmin_tmax = each_taper_tmin_tmax
            min_period, max_period = map(float, each_period.split(","))
            for each_asdf_path, each_cmt_path, each_conv_path in zip(asdf_paths_this_rank, cmt_paths_this_rank, conv_paths_this_rank):
                cmtsolution = obspy.read_events(each_cmt_path)[0]
                tau = cmtsolution.focal_mechanisms[0].moment_tensor.source_time_function.duration / 2
                input_virasdf = VirAsdf()
                input_virasdf.read_asdf(each_asdf_path)
                output_virasdf = conv_process_single_virasdf(input_virasdf, waveform_length, tau, 0,
                                                             taper_tmin_tmax, min_period, max_period, sampling_rate)
                # here each_conv_path is not the form we are using in processing sync, so do some modifications here
                each_conv_path = each_conv_path + \
                    f"preprocessed_{int(min_period)}s_to_{int(max_period)}s.h5"
                output_virasdf.write_asdf(each_conv_path, "convolve_processed")
    main()
