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
        if (asdf_path == None):
            return
        if (isfile(asdf_path)):
            # since virasdf will not used parallel io, we disable mpi here
            with pyasdf.ASDFDataSet(asdf_path, mode="r", mpi=False) as ds:
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
        # since virasdf will not used parallel io, we disable mpi here
        with pyasdf.ASDFDataSet(output_path, mode="w", mpi=False) as ds:
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
        self.events = events.copy()

    def set_waveforms(self, waveforms):
        self.waveforms = waveforms
        self.waveforms_list = list(self.waveforms.keys())

    def __copy__(self):
        new_virasdf = VirAsdf()
        new_virasdf.set_events(self.events)
        new_waveforms = {}
        for each_net_sta in self.waveforms_list:
            new_waveforms[each_net_sta] = {
                # as we don't modify inv or from the dict this value is a copy?
                "inv": self.waveforms[each_net_sta]["inv"],
                "st": self.waveforms[each_net_sta]["st"].copy()
            }
        new_virasdf.set_waveforms(new_waveforms)
        return new_virasdf

    def update_st(self, st, net_sta):
        if (self.waveforms[net_sta]["st"]):
            self.waveforms[net_sta]["st"] = st
