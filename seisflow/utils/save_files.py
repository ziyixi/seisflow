import pickle
from os.path import basename, join


def save_pickle_event(to_save, output_dir, used_gcmtid):
    output_path = join(output_dir, f"{used_gcmtid}.pkl")
    with open(output_path, "wb") as f:
        pickle.dump(to_save, f)
