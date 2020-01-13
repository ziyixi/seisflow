"""
process sync file. 
"""
from os.path import join

# fix a bug in intel
import mpi4py
import numpy as np
import obspy
import pyasdf
from mpi4py import MPI
from obspy.signal.util import _npts2nfft

from ..tasks.process.process_sync_single_st import process_sync_single_trace

mpi4py.rc.recv_mprobe = False


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
isroot = (rank == 0)


def process_single_event(min_periods, max_periods, taper_tmin_tmax, asdf_filename, waveform_length, sampling_rate, output_directory):
    # with pyasdf.ASDFDataSet(asdf_filename) as ds:
    with pyasdf.ASDFDataSet(asdf_filename, mode="r") as ds:
        # some parameters
        event = ds.events[0]
        origin = event.preferred_origin() or event.origins[0]
        event_time = origin.time

        for min_period, max_period in zip(min_periods, max_periods):
            # inv will be None if no inv provided
            def process_function(st, inv):
                # wrap process_sync_single_trace
                st = process_sync_single_trace(
                    st, event_time, waveform_length, taper_tmin_tmax, min_period, max_period, sampling_rate)

                return st
            tag_name = "preprocessed_%is_to_%is" % (
                int(min_period), int(max_period))
            tag_map = {
                "synthetic": tag_name
            }

            output_name_head = asdf_filename.split("/")[-1].split(".")[0]
            ds.process(process_function, join(
                output_directory, output_name_head+"."+tag_name + ".h5"), tag_map=tag_map)
