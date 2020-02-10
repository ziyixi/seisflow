"""
combine_asdf_files.py: Combine two adsf files.
"""
from pyasdf import ASDFDataSet


def combine_asdf(base_asdf_path, append_asdf_path, output_asdf_path):
    """
    combine_asdf: merge the waveforms in append_asdf to base_asdf, and generate a new asdf.
    """
    base_asdf = ASDFDataSet(base_asdf_path, mode="r", mpi=False)
    append_asdf = ASDFDataSet(append_asdf_path, mode="r", mpi=False)
    output_asdf = ASDFDataSet(output_asdf_path, mpi=False)
    # * add events
    events = base_asdf.events
    event = events[0]
    output_asdf.add_quakeml(events)
    # * add waveforms and stationxml
    # firstly we add base asdf
    rep_net_sta = base_asdf.waveforms.list()[0]
    tag_default = base_asdf.waveforms[rep_net_sta].get_waveform_tags()[0]
    for each_net_sta in base_asdf.waveforms.list():
        tag = base_asdf.waveforms[each_net_sta].get_waveform_tags()[0]
        assert tag == tag_default
        st = base_asdf.waveforms[each_net_sta][tag]
        inv = base_asdf.waveforms[each_net_sta]["StationXML"]
        output_asdf.add_waveforms(st, tag=tag, event_id=event)
        output_asdf.add_stationxml(inv)
    # secondly we add append asdf
    for each_net_sta in append_asdf.waveforms.list():
        tag = append_asdf.waveforms[each_net_sta].get_waveform_tags()[0]
        assert tag == tag_default
        st = append_asdf.waveforms[each_net_sta][tag]
        inv = append_asdf.waveforms[each_net_sta]["StationXML"]
        output_asdf.add_waveforms(st, tag=tag, event_id=event)
        output_asdf.add_stationxml(inv)
    del base_asdf
    del append_asdf
    del output_asdf
