"""
hinet_extracter.py: extract hinet data to sac and pz files.
"""
from glob import glob
from os.path import join, isdir, dirname

import click
import sh
from HinetPy import win32
import numpy as np


def extract_sac(data, ctable, processes):
    """
    extract the downloaded data to the sac and pz format.
    """
    # * get dirname
    thedir = dirname(data)
    # * mkdir for the sac and pz files
    sh.mkdir("-p", join(thedir, "SAC"))
    sh.mkdir("-p", join(thedir, "PZ"))
    # * extract sac
    win32.extract_sac(data, ctable, with_pz=False, outdir=join(
        thedir, "SAC"), processes=processes)
    # * extract pz
    win32.extract_pz(ctable, outdir=join(thedir, "PZ"))


@click.command()
@click.option('--data_directory', required=True, type=str, help="the data directory")
@click.option('--processes', required=True, type=int, help="the number of processes to use")
def main(data_directory, processes):
    """
    extrct the data from the downloaded database.
    """
    all_event_directories = sorted(glob(join(data_directory, "*")))
    all_event_directories = [
        item for item in all_event_directories if isdir(item)]
    sh.touch(join(data_directory, "extract.filelist"))
    extreacted_event_directories = np.loadtxt(
        join(data_directory, "extract.filelist"), dtype=np.str)
    to_extract_event_directories = sorted(
        set(all_event_directories)-set(extreacted_event_directories))
    for each_data_directory in to_extract_event_directories:
        ctable = glob(join(each_data_directory, "*ch"))[0]
        data = glob(join(each_data_directory, "*cnt"))[0]
        extract_sac(data, ctable, processes)
        with open(join(data_directory, "extract.filelist"), "a") as file:
            file.write(each_data_directory+"\n")


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
