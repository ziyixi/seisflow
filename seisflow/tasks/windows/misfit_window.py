import configparser
import sys
import warnings

import numpy as np
import pyasdf
from obspy.signal.cross_correlation import correlate, xcorr_max

from .window import Window

if not sys.warnoptions:
    warnings.simplefilter("ignore")


class Misfit_window(Window):
    def __init__(self, parent_window):
        super().__init__(left=parent_window.left, right=parent_window.right, channel=parent_window.channel,
                         network=parent_window.network, station=parent_window.station, phases=parent_window.phases.copy(), gcmtid=parent_window.gcmtid)
        self.net_sta = f"{self.network}.{self.station}"
        self.cc = None
        self.similarity = None
        self.snr_energy = None
        self.snr_amp = None
        self.deltat = None
        # we should always keep only one channel
        self.component = self.channel[-1]

    def update_first_arrival_baz(self, first_arrival_dict, baz_dict):
        # load from file, just an float
        self.first_arrival = first_arrival_dict[self.gcmtid][self.net_sta]
        self.baz = baz_dict[self.gcmtid][self.net_sta]

    def update_snr(self, data_asdf, sync_asdf):
        if (self.net_sta not in data_asdf.waveforms.list()):
            return
        data_wg = data_asdf.waveforms[self.net_sta]
        data_tag = data_wg.get_waveform_tags()[0]
        data_tr = data_wg[data_tag].select(component=self.component)[0].copy()
        event_time = sync_asdf.events[0].origins[0].time
        # get the noise window
        signal_st = data_tr.slice(self.left, self.right)
        signal_st.taper(0.05, type="hann")
        # get averaged power ratio
        signal_data = signal_st.data
        # noise_avg_power = np.sum(noise_data**2) / len(noise_data)
        noise_avg_power, noise_max_amp = self.cal_noise_average_energy(
            data_tr, self.first_arrival, event_time)
        signal_avg_power = np.sum(signal_data ** 2) / len(signal_data)
        signal_max_amp = np.max(np.abs(signal_data))
        self.snr_energy = 10 * np.log10(signal_avg_power / noise_avg_power)
        self.snr_amp = 20 * \
            np.log10(np.abs(signal_max_amp) / np.abs(noise_max_amp))

    def update_cc_deltat(self, data_asdf, sync_asdf):
        if (self.net_sta not in data_asdf.waveforms.list()):
            return
        # we assume the delta and the event_time are the same, but the starttime may be slightly different
        # also we have to make sure the net_sta is existing
        data_wg = data_asdf.waveforms[self.net_sta]
        data_tag = data_wg.get_waveform_tags()[0]
        sync_wg = sync_asdf.waveforms[self.net_sta]
        sync_tag = sync_wg.get_waveform_tags()[0]
        data_tr = data_wg[data_tag].select(component=self.component)[0].copy()
        sync_tr = sync_wg[sync_tag].select(component=self.component)[0].copy()
        # we make the starttime of sync to be the same with data
        tolerance_time = 60
        time_difference = np.abs(
            sync_tr.stats.starttime - data_tr.stats.starttime)
        if (time_difference <= data_tr.stats.delta):
            sync_tr.stats.starttime = data_tr.stats.starttime
        elif ((time_difference <= tolerance_time) and (data_tr.stats.starttime <= self.left)):
            if(sync_tr.stats.starttime < data_tr.stats.starttime):
                sync_tr.trim(data_tr.stats.starttime, sync_tr.stats.endtime)
                sync_tr.stats.starttime = data_tr.stats.starttime
            else:
                data_tr.trim(sync_tr.stats.starttime, data_tr.stats.endtime)
                sync_tr.stats.starttime = data_tr.stats.starttime
        else:
            return
        # cut to the window
        data_win_tr = data_tr.slice(self.left, self.right)
        data_win_tr.taper(0.05, type="hann")
        sync_win_tr = sync_tr.slice(self.left, self.right)
        sync_win_tr.taper(0.05, type="hann")
        # use data as the reference, calculate cc and deltat
        cc_all = correlate(data_win_tr, sync_win_tr, None,
                           demean=False, normalize="naive")
        self.similarity = cc_all[len(cc_all) // 2]
        self.deltat, self.cc = xcorr_max(cc_all, abs_max=False)
        delta = data_tr.stats.delta
        self.deltat = self.deltat * delta

    def __repr__(self):
        return f"Windows(left={self.left},right={self.right},channel={self.channel},network={self.network},gcmtid={self.gcmtid},station={self.station},phases={self.phases},snr_energy={self.snr_energy},snr_amp={self.snr_amp},deltat={self.deltat},cc={self.cc},similarity={self.similarity})"

    def cal_noise_average_energy(self, data_tr, first_arrival, event_time):
        """
        Because the data may have a pulse at the beginning, we should remove these parts
        """
        noise_start = None
        if(first_arrival == None):
            # no first arrival, we don't use that trace
            return 1e9
        noise_start = 0
        noise_win_start = event_time + noise_start
        # avoid containning the first arrival
        if(first_arrival > 30):
            noise_win_end = event_time + first_arrival - 20
        elif (first_arrival > 20):
            noise_win_end = event_time + first_arrival - 10
        else:
            noise_win_end = event_time+first_arrival
        tr_noise = data_tr.slice(noise_win_start, noise_win_end)
        # tr_noise.taper(0.05, type="hann")
        noise_max_amp = np.max(np.abs(tr_noise.data))
        noise_average_energy = np.sum(tr_noise.data**2)/len(tr_noise.data)
        return noise_average_energy, noise_max_amp
