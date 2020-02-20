"""
sac2asdf_hinet.py: convert the sac files to the asdf format, not storing the response information. (as it's PZ)
"""
from glob import glob
from os.path import basename, join

import numpy as np
import obspy
import pyasdf
from obspy.core.inventory import Channel, Inventory, Network, Station


def sac2asdf_hinet(sac_directory, cmt_path, output_path):
    with pyasdf.ASDFDataSet(output_path, mode="w", compression=None, mpi=False) as ds:
        # read in eventxml
        event_xml = obspy.read_events(cmt_path)
        # add eventxml to ds
        ds.add_quakeml(event_xml)
        event = ds.events[0]
        # read in waves
        files = sorted(glob(join(sac_directory, "*")))
        inv = Inventory()  # pylint: disable=no-value-for-parameter
        net_inv = Network(code="N", stations=[])
        # * we should sort files based on the station names
        sta_collection = {}
        # * here we add the waveforms along with the process of building the inventory
        for each_file in files:
            tr = obspy.read(each_file)[0]
            # here we need to modify some stats' values
            net_sta = tr.stats.station
            net, sta = net_sta.split(".")
            tr.stats.network = net
            tr.stats.station = sta
            # we change the channel names U->HHZ N->HHN E->HHE
            channel_mapper = {
                "U": "HHZ",
                "N": "HHN",
                "E": "HHE"
            }
            try:
                tr.stats.channel = channel_mapper[tr.stats.channel]
            except KeyError:
                continue
            # we have to consider the time difference in Japan
            tr.stats.starttime = tr.stats.starttime-9*60*60
            # * add the waveforms
            tr.data = np.require(tr.data, dtype="float32")
            ds.add_waveforms(tr, tag="raw", event_id=event)
            # * handle the stationxml
            cha = Channel(
                code=tr.stats.channel,
                location_code="",
                latitude=tr.stats.sac.stla,
                longitude=tr.stats.sac.stlo,
                elevation=tr.stats.sac.stel,
                depth=0.0,
                sample_rate=tr.stats.sampling_rate)
            if (sta in sta_collection):
                sta_collection[sta].channels.append(cha)
            else:
                sta_collection[sta] = Station(
                    code=sta,
                    latitude=tr.stats.sac.stla,
                    longitude=tr.stats.sac.stlo,
                    elevation=tr.stats.sac.stel)
                sta_collection[sta].channels.append(cha)
        # * now we can add all the sta to net
        for sta in sta_collection:
            if(len(sta_collection[sta].channels) == 3):
                net_inv.stations.append(sta_collection[sta])
        # * we can add net to station_xml
        inv.networks.append(net_inv)
        # * now we can add inv to asdf
        ds.add_stationxml(inv)
