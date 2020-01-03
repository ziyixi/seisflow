import pickle
from glob import glob
from os.path import basename, join

import numpy as np
from mpi4py import MPI
from pyasdf import ASDFDataSet

from ..asdf.extract_data_info import extract_data_info

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

phase_list = ["S", "sS", "SS", "P",
              "pP", "sP", "PP", "3.3kmps", "4.6kmps", "ScS"]
extra_save_list = ["event_time", "gcarc",
                   "az", "baz", "evla", "evlo", "evdp"]


def load_asdf_info(asdf_fname):
    # asdf file
    with ASDFDataSet(asdf_fname, mode="r") as asdf_file:
        lat = asdf_file.events[0].preferred_origin().latitude
        lon = asdf_file.events[0].preferred_origin().longitude
        dep = asdf_file.events[0].preferred_origin().depth
        time = asdf_file.events[0].preferred_origin().time

    return lat, lon, dep, time


def load_station_info(station_fname):
    # station file
    stations = np.loadtxt(station_fname, dtype=np.str)
    return stations


def get_used_asdf_files_in_directory(asdf_directory):
    all_files = glob(join(asdf_directory, "*h5"))
    used_files = np.array_split(all_files, size)[rank]
    return used_files


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--asdf_directory', required=True, type=str, help="asdf files directory.")
    @click.option('--station_fname', required=True, type=str, help="the station fname in Specfem format.")
    @click.option('--output_dir', required=True, type=str, help="the output directory to store the pickle files.")
    def main(asdf_directory, station_fname, output_dir):
        used_files = get_used_asdf_files_in_directory(asdf_directory)
        result_this_rank = {}
        for each_file in used_files:
            gcmtid = basename(each_file).split(".")[0]
            lat, lon, dep, time = load_asdf_info(each_file)
            stations = load_station_info(station_fname)
            result_this_rank[gcmtid] = extract_data_info(
                gcmtid, lat, lon, dep, time, stations)
        comm.Barrier()
        result_gathered = comm.gather(result_this_rank, root=0)
        # collect all results into different pickle files
        if(rank == 0):
            for each_phase in phase_list:
                output_result = {}
                output_fname = join(output_dir, f"traveltime.{each_phase}.pkl")
                for result_each_rank in result_gathered:
                    for gcmtid in result_each_rank:
                        output_result[gcmtid] = {}
                        for net_sta in result_each_rank[gcmtid]:
                            travel_time = result_each_rank[gcmtid][net_sta][each_phase]
                            output_result[gcmtid][net_sta] = travel_time
                with open(output_fname, "wb") as handle:
                    pickle.dump(output_result, handle)
            for each_extra in extra_save_list:
                output_result = {}
                output_fname = join(output_dir, f"extra.{each_extra}.pkl")
                for result_each_rank in result_gathered:
                    for gcmtid in result_each_rank:
                        output_result[gcmtid] = {}
                        for net_sta in result_each_rank[gcmtid]:
                            extra_value = result_each_rank[gcmtid][net_sta][each_extra]
                            output_result[gcmtid][net_sta] = extra_value
                with open(output_fname, "wb") as handle:
                    pickle.dump(output_result, handle)
    main()
