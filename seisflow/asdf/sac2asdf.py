"""
sac2asdf: convert sac files to the asdf format. Station Response file should be in RESP format.
"""
import subprocess
import tempfile
from glob import glob
from os.path import join

import obspy
import pyasdf


def sac2asdf(sac_directory, response_directory, cmt_path, output_path):
    with pyasdf.ASDFDataSet(output_path, mode="w", compression=None) as ds:
        # read in eventxml
        event_xml = obspy.read_events(cmt_path)
        # add eventxml to ds
        ds.add_quakeml(event_xml)
        event = ds.events[0]
        # read in waves
        files = glob(join(sac_directory, "*"))
        station_xml = obspy.core.inventory.inventory.Inventory()

        for i, filename in enumerate(files):
            waveform_stream = obspy.read(filename)
            ds.add_waveforms(waveform_stream, tag="raw", event_id=event)

            # add stationxml
            allfiles = glob(join(response_directory, "*"))
            for fname in allfiles:
                station_xml_this_seed = obspy.core.inventory.inventory.Inventory()
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
