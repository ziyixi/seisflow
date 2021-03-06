"""
hinet_downloader.py: download the hinet data based on HinetPy.
"""
from datetime import timedelta
from glob import glob
from os.path import basename, dirname, join

import click
import numpy as np
import obspy
import sh
from HinetPy import Client, win32
from loguru import logger


def init_client(username, password):
    """
    init client and check status.
    """
    client = Client(username, password)
    client.doctor()
    return client


def download_file(each_cmt_file, client, output_directory, network_code):
    """
    download for a single event.
    """
    # * mkdir for the event
    cmt = basename(each_cmt_file)
    sh.mkdir("-p", join(output_directory, cmt))
    # * load the cmt info
    event = obspy.read_events(each_cmt_file)[0]
    origin = event.preferred_origin() or event.origins[0]
    eventtime = origin.time.datetime + timedelta(hours=9)
    starttime = eventtime - timedelta(minutes=2)
    try:
        data, ctable = client.get_continuous_waveform(
            network_code, starttime, 42, outdir=join(output_directory, cmt))
    except ValueError:
        return None, None
    return data, ctable


@click.command()
@click.option('--logfile', required=True, type=str, help="the log file path")
@click.option('--cmts_directory', required=True, type=str, help="the cmts directory path")
@click.option('--output_directory', required=True, type=str, help="the output directory path")
@click.option('--username', required=True, type=str, help="the username for the hinet website.")
@click.option('--password', required=True, type=str, help="the password for the hinet website.")
@click.option('--network_code', required=True, type=str, help="the network code, eg: 0101 for hi-net and 0120 for s-net")
def main(logfile, cmts_directory, output_directory, username, password, network_code):
    """
    Download Hinet waveforms based on cmtsolution files.
    """
    # * set up logging
    logger.add(logfile, format="{time} {level} {message}", level="INFO")
    # * start the client
    client = init_client(username, password)
    logger.info("start to download.")
    # * here we get the cmtsolution files paths
    all_cmt_files = sorted(glob(join(cmts_directory, "*")))
    # * mkdir for the data
    sh.mkdir("-p", output_directory)
    # * in the output directory, we use a file called local.filelist to select the cmt files to download
    sh.touch(join(output_directory, "local.filelist"))
    installed_cmt_files = np.loadtxt(
        join(output_directory, "local.filelist"), dtype=np.str)
    to_download_cmt_files = sorted(
        set(all_cmt_files) - set(installed_cmt_files))
    data_list = []
    ctable_list = []
    for each_cmt_file in to_download_cmt_files:
        data, ctable = download_file(
            each_cmt_file, client, output_directory, network_code)
        if((data != None) and (ctable != None)):
            data_list.append(data)
            ctable_list.append(ctable)
            with open(join(output_directory, "local.filelist"), 'a') as file:
                file.write(each_cmt_file + "\n")
            logger.success(f"download {basename(each_cmt_file)}")
        else:
            logger.error(f"skip {basename(each_cmt_file)}")


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
