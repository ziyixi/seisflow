"""
generate_adjoint_source_for_single_event_station_pair: generate the adjoint source for se, used to check the kernel in the structure inversion.
"""
from ...utils.asdf_io import VirAsdf
from .adjoint_source_each_window_zerolagcc import \
    calculate_adjoint_source_each_window
import numpy as np
from ...utils.load_files import load_pickle


def generate_adjoint_source_trace(net_sta, category, phases, misfit_windows_all, raw_sync_virasdf, sync_viradsf, data_virasdf, mintime, maxtime):
    """
    Prepare the parameters used to generate the adjoint source for the event-station pair.
    """
    # get mistit_window
    misfit_windows_net_sta_category = misfit_windows_all[net_sta][category].windows
    # * according to the phase, select the misfit window
    # to prevent the possible problem, the phases is required to be the same with the window
    misfit_windows_used = None
    for each_misfit_window in misfit_windows_net_sta_category:
        if (phases in each_misfit_window.phases):
            misfit_windows_used = each_misfit_window
    if (misfit_windows_used == None):
        raise Exception("check the phases provided!")
    print("="*100)
    print("information about the window used:")
    print(misfit_windows_used)
    print("="*100)
    # * get raw_sync_asdf_trace
    # here we only consider the one window case.
    component = misfit_windows_used.channel[-1]
    st_raw = raw_sync_virasdf.get_waveforms()[net_sta]["st"]
    raw_sync_asdf_trace_used = st_raw.select(component=component)[0]
    # * get sync_asdf_trace
    st_sync = sync_viradsf.get_waveforms()[net_sta]["st"]
    sync_asdf_trace_used = st_sync.select(component=component)[0]
    # * get data_asdf_trace
    st_data = data_virasdf.get_waveforms()[net_sta]["st"]
    data_asdf_trace_used = st_data.select(component=component)[0]
    adjoint_source_trace = calculate_adjoint_source_each_window(
        misfit_windows_used, raw_sync_asdf_trace_used, sync_asdf_trace_used, data_asdf_trace_used, mintime, maxtime)
    return adjoint_source_trace, misfit_windows_used


def convert_adjoint_source_to_array(misfit_windows_used, adjoint_source_trace, category):
    """
    convert the adjoint source in the trace format to the array required by the adjoint source.
    """
    adjoint_source_length = len(adjoint_source_trace.data)
    adjoint_source_array = np.zeros((3, adjoint_source_length))
    component = category[-1]
    if (component == "z"):
        adjoint_source_array[2, :] = adjoint_source_trace.data[:]
    elif (component == "r"):
        theta = (misfit_windows_used.baz - 180) % 360
        adjoint_source_array[0, :] = adjoint_source_trace.data[:] * \
            np.sin(np.deg2rad(theta))
        adjoint_source_array[1, :] = adjoint_source_trace.data[:] * \
            np.cos(np.deg2rad(theta))
    elif (component == "t"):
        theta = (misfit_windows_used.baz - 90) % 360
        adjoint_source_array[0, :] = adjoint_source_trace.data[:] * \
            np.sin(np.deg2rad(theta))
        adjoint_source_array[1, :] = adjoint_source_trace.data[:] * \
            np.cos(np.deg2rad(theta))
    else:
        raise Exception("no such component!")
    return adjoint_source_array


if __name__ == "__main__":
    import click
    import pyasdf

    def save_adjoint_to_asdf(net_sta, adjoint_source_array, output_fname):
        """
        save_adjoint_to_asdf: save the adjoint source to asdf format to be read by Specfem.
        """
        components = ["MXE", "MXN", "MXZ"]
        with pyasdf.ASDFDataSet(output_fname, mode="w", mpi=False) as output_asdf:
            for index_component in range(3):
                component = components[index_component]
                specfem_adj_source = adjoint_source_array[index_component, :]
                tag = net_sta.replace(".", "_") + "_" + \
                    components[index_component]
                output_asdf.add_auxiliary_data(
                    data=specfem_adj_source, data_type="AdjointSources", path=tag, parameters={})

    @click.command()
    @click.option('--net_sta', required=True, type=str, help="network and station name, eg: net.sta")
    @click.option('--category', required=True, type=str, help="the category used, can be z, r, t, surface_z, surface_r and surface_t")
    @click.option('--phases', required=True, type=str, help="the phases of the window, eg: P")
    @click.option('--filter_time_range', required=True, type=str, help="the time range to do the bandpass filter, eg: 10,40")
    @click.option('--misfit_windows_path', required=True, type=str, help="the path of the misfit windows pickle file")
    @click.option('--raw_sync_path', required=True, type=str, help="the path of the raw sync asdf")
    @click.option('--sync_path', required=True, type=str, help="the path of the processed sync asdf")
    @click.option('--data_path', required=True, type=str, help="the path of the processed data asdf")
    @click.option('--output_path', required=True, type=str, help="the output path of the adjoint source")
    def main(net_sta, category, phases, filter_time_range, misfit_windows_path, raw_sync_path, sync_path, data_path, output_path):
        # * firstly we load the asdf files to Virasdf
        raw_sync_virasdf = VirAsdf()
        raw_sync_virasdf.read_asdf(raw_sync_path)
        sync_viradsf = VirAsdf()
        sync_viradsf.read_asdf(sync_path)
        data_virasdf = VirAsdf()
        data_virasdf.read_asdf(data_path)
        # * we can generate the adjoint source trace
        misfit_windows_all = load_pickle(misfit_windows_path)
        mintime, maxtime = map(float, filter_time_range.split(","))
        adjoint_source_trace, misfit_windows_used = generate_adjoint_source_trace(
            net_sta, category, phases, misfit_windows_all, raw_sync_virasdf, sync_viradsf, data_virasdf, mintime, maxtime)
        # * convert the trace to the numpy array
        adjoint_source_array = convert_adjoint_source_to_array(
            misfit_windows_used, adjoint_source_trace, category)
        # * save the adjoint source
        gcmtid = misfit_windows_used.gcmtid
        save_adjoint_to_asdf(net_sta, adjoint_source_array,
                             output_path)

    main()
