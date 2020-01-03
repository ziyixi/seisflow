import subprocess

import click
import obspy
import pyasdf


def remove_stations(ds):
    waveform_list = ds.waveforms.list()
    for item in waveform_list:
        wg = ds.waveforms[item]
        data_tags = wg.get_waveform_tags()
        if(len(data_tags) == 0):
            inv = wg["StationXML"]
            del ds.waveforms[item]
        else:
            # see if there are only three traces
            tag = data_tags[0]
            st = wg[tag]
            # after procesing the data, the only possible is the case of mul locs
            if(len(st) == 3):
                continue
            else:
                loc_set = set()
                for trace in st:
                    net, sta, loc, cha = trace.id.split(".")
                    loc_set.add(loc)
                if("" in loc_set):
                    reference_loc = ""
                else:
                    reference_loc = sorted(loc_set)[0]

                for trace in st:
                    net, sta, loc, cha = trace.id.split(".")
                    if(loc != reference_loc):
                        st.remove(trace)

                inv = wg["StationXML"]
                discard_locs = sorted(loc_set-set([reference_loc]))


def simplify_data_asdf(asdf_file):
    with pyasdf.ASDFDataSet(asdf_file) as ds:
        remove_stations(ds)
