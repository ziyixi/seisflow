import pickle
from glob import glob
from os.path import basename, join

import click
import numpy as np
import tqdm

from ...tasks.adjoint.calculate_adjoint_source_zerolag_cc_multiple_events import (
    get_weights_for_all, update_categorical_weighting_multiple_events)
from ...utils.load_files import load_pickle


@click.command()
@click.option('--misfit_directory', required=True, type=str, help="the misfit windows directory")
@click.option('--output_path', required=True, type=str, help="the output weight path")
@click.option('--stations_path', required=True, type=str, help="the stations path")
@click.option('--snr', required=True, type=str, help="snr1,snr2")
@click.option('--cc', required=True, type=str, help="cc1,cc2")
@click.option('--deltat', required=True, type=str, help="deltat1,deltat2")
def main(misfit_directory, output_path, stations_path, snr, cc, deltat):
    all_paths = sorted(glob(join(misfit_directory, "*pkl")))
    misfit_windows_dict = {}
    for each_path in all_paths:
        gcmtid = basename(each_path).split(".")[0]
        misfit_windows_dict[gcmtid] = load_pickle(each_path)
    stations = np.loadtxt(stations_path, dtype=np.str)
    snr_threshold = tuple(map(float, snr.split(",")))
    cc_threshold = tuple(map(float, cc.split(",")))
    deltat_threshold = tuple(map(float, deltat.split(",")))

    weight_dict = {}
    for gcmtid in tqdm.tqdm(misfit_windows_dict):
        misfit_windows = misfit_windows_dict[gcmtid]
        weight_dict[gcmtid] = get_weights_for_all(
            misfit_windows, stations,  snr_threshold, cc_threshold, deltat_threshold, False, print_info=False)
    weight_dict = update_categorical_weighting_multiple_events(weight_dict)
    with open(output_path, "wb") as f:
        pickle.dump(weight_dict, f)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
