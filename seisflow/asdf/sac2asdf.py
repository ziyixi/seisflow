"""
sac2asdf: convert sac files to the asdf format. 
"""
from glob import glob
from os.path import join, isfile

import obspy
import pyasdf
from obspy.core.inventory.inventory import Inventory


def sac2asdf(sac_directory, response_directory, cmt_path, output_path):
    with pyasdf.ASDFDataSet(output_path, mode="w", compression=None, mpi=False) as ds:
        # read in eventxml
        event_xml = obspy.read_events(cmt_path)
        # add eventxml to ds
        ds.add_quakeml(event_xml)
        event = ds.events[0]
        # read in waves
        files = sorted(glob(join(sac_directory, "*")))
        # we should select files
        files = [item for item in files if isfile(item)]
        station_xml = Inventory()  # pylint: disable=no-value-for-parameter
        stream_mapper = {}
        for filename in files:
            waveform_stream = obspy.read(filename)
            waveid = waveform_stream[0].id
            net, sta, _, _ = waveid.split(".")
            thekey = f"{net}.{sta}"
            stream_mapper[thekey] = waveform_stream
            ds.add_waveforms(waveform_stream, tag="raw", event_id=event)

        #     # add stationxml
        #     allfiles = sorted(glob(join(response_directory, "*")))
        #     station_xml_this_seed = Inventory()  # pylint: disable=no-value-for-parameter
        #     for fname in allfiles:
        #         inv_temp = obspy.read_inventory(fname)
        #         # update essencial location information
        #         inv_temp = update_info(inv_temp, waveform_stream)
        #         if(inv_temp == None):
        #             continue
        #         station_xml_this_seed += inv_temp

        #     station_xml += station_xml_this_seed

        # ds.add_stationxml(station_xml)

        # add stationxml
        allfiles = sorted(glob(join(response_directory, "*")))
        for fname in allfiles:
            inv_temp = obspy.read_inventory(fname)
            rep_usable_channel = inv_temp.get_contents()["channels"][0]
            net, sta, _, _ = rep_usable_channel.split(".")
            thekey = f"{net}.{sta}"
            waveform_stream = stream_mapper[thekey]
            inv_temp = update_info(inv_temp, waveform_stream)
            station_xml += inv_temp
        ds.add_stationxml(station_xml)


def update_info(inv, waveform_stream):
    usable_channels = inv.get_contents()["channels"]
    # loop all channels, search info
    status = False
    for channel in usable_channels:  # channel, like BO.NKG..BHE
        if(status == False):
            for thewave in waveform_stream:
                waveid = thewave.id
                waveid_split = waveid.split(".")
                channel_split = channel.split(".")
                waveid_key = ".".join(
                    [waveid_split[0], waveid_split[1]])
                channel_key = ".".join(
                    [channel_split[0], channel_split[1]])
                # we assume the location id will not change the lle
                if(waveid_key == channel_key):
                    status = True
                    # note here we only update the station part of the inv, use with cautious!
                    inv[0][0].latitude = thewave.stats.sac.stla
                    inv[0][0].longitude = thewave.stats.sac.stlo
                    inv[0][0].elevation = thewave.stats.sac.stel
                    break
    if(not status):
        return None
    else:
        return inv
