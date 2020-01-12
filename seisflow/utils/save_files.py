import pickle
from os.path import basename, join


def save_pickle_event(to_save, output_dir, used_gcmtid):
    output_path = join(output_dir, f"{used_gcmtid}.pkl")
    with open(output_path, "wb") as f:
        pickle.dump(to_save, f)


def save_cmtsolution(output_path, cmtsolution):
    """
    save cmtsolution file.
    """
    cmtsolution.write(output_path, format="CMTSOLUTION")
    # remove the final line (or specfem will report error)
    with open(output_path, "r") as f:
        lines = f.readlines()
    with open(output_path, "w") as f:
        for i in range(12):
            f.write(lines[i])
        f.write(lines[-1].split("\n")[0])
