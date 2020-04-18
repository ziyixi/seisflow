"""
tao_2018_ggg_windows: get windows described in tao et al. ggg, 2018. Only works for one event (gcmtid).
"""
from os.path import join

import numpy as np
from mpi4py import MPI

from ...tasks.windows.tao_2018_ggg_windows import generate_windows
from ...utils.load_files import load_pickle
from ...utils.save_files import save_pickle_event

phase_list = ["S", "sS", "SS", "P",
              "pP", "sP", "PP", "3.3kmps", "4.6kmps", "ScS"]
comm = MPI.COMM_WORLD  # pylint: disable=c-extension-no-member
size = comm.Get_size()
rank = comm.Get_rank()
# * test is passed for this script on 01/07/2020


def load_traveltime(data_info_directory):
    """
    load the traveltime pickle files and convert it to traveltime_all[gcmtid][net_sta][each_phase]
    """
    # load files
    loaded_dicts = {}
    for each_phase in phase_list:
        fname = join(data_info_directory, f"traveltime.{each_phase}.pkl")
        loaded_dicts[each_phase] = load_pickle(fname)
    # convert to traveltime_all[gcmtid][net_sta][each_phase]
    traveltime_all = {}
    all_gcmtids = list(loaded_dicts[phase_list[0]].keys())
    all_net_stas = list(loaded_dicts[phase_list[0]][all_gcmtids[0]].keys())
    for each_gcmtid in all_gcmtids:
        traveltime_all[each_gcmtid] = {}
        for each_net_sta in all_net_stas:
            traveltime_all[each_gcmtid][each_net_sta] = {}
            for each_phase in phase_list:
                traveltime_all[each_gcmtid][each_net_sta][each_phase] = loaded_dicts[each_phase][each_gcmtid][each_net_sta]
    return traveltime_all


def load_eventtime(data_info_directory):
    fpath = join(data_info_directory, f"extra.event_time.pkl")
    eventtime = load_pickle(fpath)
    return eventtime


def get_traveltime_all_this_rank(traveltime_all):
    traveltime_all_this_rank = {}
    all_gcmtids = sorted(list(traveltime_all.keys()))
    all_gcmtids_this_rank = np.array_split(all_gcmtids, size)[rank]
    for each_gcmtid_this_rank in all_gcmtids_this_rank:
        traveltime_all_this_rank[each_gcmtid_this_rank] = traveltime_all[each_gcmtid_this_rank]
    return traveltime_all_this_rank


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--data_info_directory', required=True, type=str, help="the data info directory from extract data info script")
    @click.option('--time_length', required=True, type=float, help="maximum time length for the windows")
    @click.option('--output_dir', required=True, type=str, help="the output directory to store the window files")
    def main(data_info_directory, time_length, output_dir):
        traveltime_all = load_traveltime(data_info_directory)
        eventtime = load_eventtime(data_info_directory)
        traveltime_all_this_rank = get_traveltime_all_this_rank(traveltime_all)
        for each_gcmtid in traveltime_all_this_rank:
            each_traveltime = traveltime_all_this_rank[each_gcmtid]
            rep_net_sta = list(eventtime[each_gcmtid].keys())[0]
            each_event_time = eventtime[each_gcmtid][rep_net_sta]
            windows_this_gcmtid = generate_windows(
                each_gcmtid, each_traveltime, each_event_time, time_length)
            save_pickle_event(windows_this_gcmtid, output_dir, each_gcmtid)

    main()  # pylint: disable=no-value-for-parameter
