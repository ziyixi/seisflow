import pickle
from os.path import join, basename
from .load_files import load_pickle


def save_pickle_event(to_save, output_dir, used_gcmtid):
    output_path = join(output_dir, f"{used_gcmtid}.pkl")
    with open(output_path, "wb") as f:
        pickle.dump(to_save, f)


def load_windows(gcmtid, windows_directory):
    fname = join(windows_directory, f"{gcmtid}.pkl")
    windows = load_pickle(fname)
    return windows
