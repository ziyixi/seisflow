import pickle
from os.path import basename, join


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


def load_gcarc(data_info_directory):
    gcarc_path = join(
        data_info_directory, "extra.gcarc.pkl")
    gcarc = load_pickle(gcarc_path)
    return gcarc


def get_stations_mapper(stations):
    stations_mapper = {}
    for row in stations:
        net = row[1]
        sta = row[0]
        net_sta = f"{net}.{sta}"
        # lon,lat
        stations_mapper[net_sta] = (float(row[3]), float(row[2]))
    return stations_mapper


def load_pickle_event(pickle_directory, used_gcmtid):
    return load_pickle(join(pickle_directory, f"{used_gcmtid}.pkl"))
