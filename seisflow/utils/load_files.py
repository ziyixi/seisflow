import pickle
from os.path import join, basename


def load_pickle(pickle_path):
    with open(pickle_path, "rb") as f:
        data = pickle.load(f)
    return data


def load_windows(gcmtid, windows_directory):
    fname = join(windows_directory, f"{gcmtid}.pkl")
    windows = load_pickle(fname)
    return windows


def load_first_arrival_baz_evdp(data_info_directory):
    first_arrival_zr_path = join(
        data_info_directory, "traveltime.P.pkl")
    first_arrival_t_path = join(
        data_info_directory, "traveltime.S.pkl")
    baz_path = join(
        data_info_directory, "extra.baz.pkl")
    evdp_path = join(
        data_info_directory, "extra.evdp.pkl")
    first_arrival_zr = load_pickle(first_arrival_zr_path)
    first_arrival_t = load_pickle(first_arrival_t_path)
    baz = load_pickle(baz_path)
    evdp = load_pickle(evdp_path)
    return first_arrival_zr, first_arrival_t, baz, evdp
