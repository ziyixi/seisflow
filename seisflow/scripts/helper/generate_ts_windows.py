"""
generate_ts_windows.py: read in the current windows pickle file and extract only the windows for the SH wave in the tangential component.
"""
from copy import deepcopy

import click

from ...utils.load_files import load_pickle
from ...utils.save_files import save_pickle


@click.command()
@click.option('--windows_path', required=True, type=str, help="the raw windows path")
@click.option('--output_path', required=True, type=str, help="the output windows path")
def main(windows_path, output_path):
    windows_raw = load_pickle(windows_path)
    windows_new = {}
    for net_sta in windows_raw:
        windows_new[net_sta] = {}
        for each_category in windows_raw[net_sta]:
            if (each_category != "t"):
                windows_new[net_sta][each_category] = windows_raw[net_sta][each_category]
                windows_new[net_sta][each_category].windows = []
            else:
                thewindow = None
                for each_window in windows_raw[net_sta][each_category].windows:
                    if (("s" in each_window.phases) or ("S" in each_window.phases)):
                        thewindow = deepcopy(each_window)
                        break
                windows_new[net_sta][each_category] = windows_raw[net_sta][each_category]
                windows_new[net_sta][each_category].windows = [thewindow]
    save_pickle(windows_new, output_path)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
