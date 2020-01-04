"""
get_stations_adjoint.py: get STATIONS_ADJOINT file used for the kernel simulation.
"""
from glob import glob
from os.path import basename, join

import numpy as np

from ..utils.load_files import load_pickle


def load_stations(stations_path):
    return np.loadtxt(stations_path, dtype=np.str)


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--stations_path', required=True, type=str, help="stations path in Specfem format")
    @click.option('--misfit_windows_directory', required=True, type=str, help="misfit windows directory")
    @click.option('--output_directory', required=True, type=str, help="output STATIONS_ADJOINT file's path")
    def main(stations_path, misfit_windows_directory, output_directory):
        stations = load_stations(stations_path)
        all_files = glob(join(misfit_windows_directory, "*pkl"))
        for each_file in all_files:
            gcmtid = basename(each_file).split(".")[0]
            each_misfit_windows = load_pickle(each_file)
            used_net_sta = list(each_misfit_windows.keys())
            output_path = join(output_directory, f"{gcmtid}")
            with open(output_path, "w") as f:
                for row in stations:
                    net_sta = f"{row[1]}.{row[0]}"
                    if (net_sta in used_net_sta):
                        f.write(" ".join(row)+"\n")
    main()
