"""
check_stations_adjoint_adjoint_source.py: check if the stations_adjoint file is compatiable with the adjoint source asdf file.
"""
from glob import glob
from os.path import join

import click
import numpy as np
import pyasdf
import tqdm


def kernel(stations_adjoint_fname, asdf_fname):
    asdf = pyasdf.ASDFDataSet(asdf_fname, mode="r")
    asdf_list = asdf.auxiliary_data.AdjointSources.list()
    stations_adjoint = np.loadtxt(stations_adjoint_fname, dtype=np.str)
    for each_asdf_net_sta in asdf_list:
        net, sta, _ = each_asdf_net_sta.split("_")
        status = False
        for row in stations_adjoint:
            if (row[0] == sta and row[1] == net):
                status = True
        if (not status):
            print(asdf_fname.split("/")[-3], net, sta)


@click.command()
@click.option('--base_directory', required=True, type=str, help="the base simulation directory")
def main(base_directory):
    all_events = sorted(glob(join(base_directory, "*")))
    for each_event in tqdm.tqdm(all_events):
        stations_adjoint_fname = join(each_event, "DATA", "STATIONS_ADJOINT")
        asdf_fname = join(each_event, "SEM", "adjoint.h5")
        kernel(stations_adjoint_fname, asdf_fname)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
