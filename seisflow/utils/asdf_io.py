"""
asdf_io.py: define the structure of virtual asdf in memory and according IO with asdf.
"""
from os.path import isfile

import numpy as np
import pyasdf


class VirAsdf():
    def __init__(self):
        super().__init__()
        self.events = None
        self.waveforms_list = []
        self.waveforms = {}
        self.asdf_path = None
        self.usecopy = False
        self.lazy = False
        self.tag = None

    def read_asdf(self, asdf_path, usecopy=False, lazy=False):
        """
        read_asdf: read in asdf
        """
        self.asdf_path = asdf_path
        self.usecopy = usecopy
        self.lazy = lazy
        if(not self.lazy):
            self.read_asdf_kernel()

    def read_asdf_kernel(self):
        if (self.asdf_path == None):
            return
        if (isfile(self.asdf_path)):
            # since virasdf will not used parallel io, we disable mpi here
            with pyasdf.ASDFDataSet(self.asdf_path, mode="r", mpi=False) as ds:
                self.events = ds.events
                self.waveforms_list = ds.waveforms.list()
                for each_net_sta in self.waveforms_list:
                    wg = ds.waveforms[each_net_sta]
                    self.tag = wg.get_waveform_tags()[0]
                    if(not self.usecopy):
                        self.waveforms[each_net_sta] = {
                            "inv": wg["StationXML"],
                            "st": wg[self.tag]
                        }
                    else:
                        self.waveforms[each_net_sta] = {
                            "inv": wg["StationXML"].copy(),
                            "st": wg[self.tag].copy()
                        }
                    # it's better to use float32
                    for tr in self.waveforms[each_net_sta]["st"]:
                        tr.data = np.require(tr.data, dtype="float32")
        else:
            pass

    def revoke(self):
        """
        revoke: revoke from the lazy status.
        """
        if (self.lazy):
            self.read_asdf_kernel()
            self.lazy = False

    def sleep(self):
        """
        sleep: remove the data that has been read.
        """
        if (not self.lazy):
            self.waveforms = None
            self.lazy = True

    def write_asdf(self, output_path, tag=None):
        """
        write_asdf: write to asdf
        """
        if (tag == None):
            tag = self.tag
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
                "inv": self.waveforms[each_net_sta]["inv"].copy(),
                "st": self.waveforms[each_net_sta]["st"].copy()
            }
        new_virasdf.set_waveforms(new_waveforms)
        new_virasdf.tag = self.tag
        return new_virasdf

    def update_st(self, st, net_sta):
        if (self.waveforms[net_sta]["st"]):
            self.waveforms[net_sta]["st"] = st
