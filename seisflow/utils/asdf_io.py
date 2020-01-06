"""
asdf_io.py: define the structure of virtual asdf in memory and according IO with asdf.
"""
from os.path import isfile

import pyasdf


class VirAsdf():
    def __init__(self):
        super().__init__()
        self.events = None
        self.waveforms_list = []
        self.waveforms = {}

    def read_asdf(self, asdf_path):
        """
        read_asdf: read in asdf
        """
        if(isfile(asdf_path)):
            with pyasdf.ASDFDataSet(asdf_path, mode="r") as ds:
                self.events = ds.events
                self.waveforms_list = ds.waveforms.list()
                for each_net_sta in self.waveforms_list:
                    wg = ds.waveforms[each_net_sta]
                    tag = wg.get_waveform_tags()[0]
                    self.waveforms[each_net_sta] = {
                        "inv": wg["StationXML"],
                        "st": wg[tag]
                    }
        else:
            pass

    def write_asdf(self, output_path, tag):
        """
        write_asdf: write to asdf
        """
        with pyasdf.ASDFDataSet(output_path, mode="w") as ds:
            ds.add_quakeml(self.events)
            for net_sta in self.waveforms:
                ds.add_waveforms(
                    self.waveforms[net_sta]["st"], tag=tag, event_id=self.events[0])
                ds.add_stationxml(self.waveforms[net_sta]["inv"])

    def get_events(self):
        return self.events

    def get_waveforms(self):
        return self.waveforms

    def get_waveforms_list(self):
        return self.waveforms_list

    def set_events(self, events):
        self.events = events

    def set_waveforms(self, waveforms):
        self.waveforms = waveforms
        self.waveforms_list = list(self.waveforms.keys())
