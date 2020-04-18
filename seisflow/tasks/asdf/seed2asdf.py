"""
seed2asdf.py: convert a seed file to an asdf file. rdseed should be installed.
"""
import subprocess
import tempfile
from glob import glob
from os.path import join

import obspy
import pyasdf
from obspy.io.xseed import Parser
from obspy.core.inventory.inventory import Inventory


def seed2asdf(seed_directory, cmt_path, output_path):
    with pyasdf.ASDFDataSet(output_path, mode="w", compression=None) as ds:
        # read in eventxml
        event_xml = obspy.read_events(cmt_path)
        # add eventxml to ds
        ds.add_quakeml(event_xml)
        event = ds.events[0]
        # read in waves
        files = sorted(glob(join(seed_directory, "*")))
        station_xml = Inventory()  # pylint: disable=no-value-for-parameter
        waveform_read_status = None

        for filename in files:
            # try to use obspy
            try:
                waveform_stream = obspy.read(filename)
                waveform_read_status = 1
            except:  # pylint: disable=bare-except
                dirpath = tempfile.mkdtemp()
                command = f"rdseed -d -f {filename} -q {dirpath}"
                subprocess.call(command, shell=True)
                waveform_stream = obspy.read(join(dirpath, "*SAC"))
                waveform_read_status = 2

            ds.add_waveforms(waveform_stream, tag="raw", event_id=event)

            # add stationxml (since statinxml may be different for different events, it's better
            # to store only one event in ds)
            try:
                station_xml_file_obj = tempfile.NamedTemporaryFile(
                    delete=False)
                station_xml_file_obj.close()
                station_xml_file_path = station_xml_file_obj.name
                sp = Parser(filename)
                sp.write_xseed(station_xml_file_path)
                station_xml_this_seed = obspy.read_inventory(
                    station_xml_file_path)
            except:  # pylint: disable=bare-except
                # since such an error occurs, we might have to use the except part to read the waveform to get the sac head
                if(waveform_read_status == 1):
                    # re readin waveform_stream
                    dirpath = tempfile.mkdtemp()
                    command = f"rdseed -d -f {filename} -q {dirpath}"
                    subprocess.call(command, shell=True)
                    waveform_stream = obspy.read(join(dirpath, "*SAC"))
                else:
                    pass  # should already have the head information

                dirpath = tempfile.mkdtemp()
                command = f"rdseed -R -f {filename} -q {dirpath}"
                subprocess.call(command, shell=True)

                station_xml_this_seed = Inventory()  # pylint: disable=no-value-for-parameter
                allfiles = sorted(glob(join(dirpath, "*")))
                for fname in allfiles:
                    inv_temp = obspy.read_inventory(fname)
                    # update essencial location information
                    inv_temp = update_info(inv_temp, waveform_stream)
                    if(inv_temp == None):
                        continue
                    station_xml_this_seed += inv_temp

            station_xml += station_xml_this_seed

        ds.add_stationxml(station_xml)


def update_info(inv, waveform_stream):
    usable_channels = inv.get_contents()["channels"]
    # loop all channels, search info
    status = False
    for channel in usable_channels:  # channel, like BO.NKG..BHE
        if(status == False):
            for thewave in waveform_stream:
                waveid = thewave.id
                if(waveid == channel):
                    status = True
                    inv[0][0].latitude = thewave.stats.sac.stla
                    inv[0][0].longitude = thewave.stats.sac.stlo
                    inv[0][0].elevation = thewave.stats.sac.stel
                    break
    if(not status):
        return None
    else:
        return inv
