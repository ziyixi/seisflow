"""
hinet_extracter.py: extract hinet data to sac and pz files.
"""
import tarfile
from glob import glob
from os.path import basename, dirname, isdir, join

import click
import numpy as np
import sh
from HinetPy import win32


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
    # * return the sac and pz directory path
    return join(thedir, "SAC"), join(thedir, "PZ")


def make_tarfile(output_filename, source_dir):
    with tarfile.open(output_filename, "w:gz") as tar:
        tar.add(source_dir, arcname=basename(source_dir))


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
        ctables = glob(join(each_data_directory, "*ch"))
        datas = glob(join(each_data_directory, "*cnt"))
        if (len(ctables) == 0 or len(datas) == 0):
            continue
        ctable = ctables[0]
        data = datas[0]
        sac_path, pz_path = extract_sac(data, ctable, processes)
        # * to save files quota, we need to zip them and remove the sac and pz files.
        # tar sac
        event_path = dirname(sac_path)
        output_sac = join(event_path, "SAC.tar.gz")
        make_tarfile(output_sac, sac_path)
        output_pz = join(event_path, "PZ.tar.gz")
        make_tarfile(output_pz, pz_path)
        # remove sac and pz
        sh.rm("-rf", sac_path)
        sh.rm("-rf", pz_path)
        with open(join(data_directory, "extract.filelist"), "a") as file:
            file.write(each_data_directory+"\n")


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
