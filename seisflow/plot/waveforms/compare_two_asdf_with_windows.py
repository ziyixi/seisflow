"""
compare_two_asdf_with_windows.py: compare two asdf files with windows.
"""
import pickle
from os.path import basename, join

import click
import matplotlib as mpl
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import numpy as np
import obspy
import pyasdf
import tqdm
from recordtype import recordtype

from ...utils.load_files import load_pickle
from ...utils.setting import SURFACE_THRESHOLD

label_size = 25
mpl.rcParams['xtick.labelsize'] = label_size

to_plot_trace = recordtype("to_plot_trace", [
                           "obs_z", "syn_z", "obs_r", "syn_r", "obs_t", "syn_t", "info"])


def getinfo(info_dir, gcmtid):
    # we load gcarc and azimuth
    gcarc_path = join(info_dir, "extra.gcarc.pkl")
    azimuth_path = join(info_dir, "extra.az.pkl")
    evdp_path = join(info_dir, "extra.evdp.pkl")
    # we load travel times
    p_path = join(info_dir, "traveltime.P.pkl")
    pp_path = join(info_dir, "traveltime.pP.pkl")
    sp_path = join(info_dir, "traveltime.sP.pkl")
    s_path = join(info_dir, "traveltime.S.pkl")
    ss_path = join(info_dir, "traveltime.sS.pkl")
    scs_path = join(info_dir, "traveltime.ScS.pkl")
    result_readin = {
        "gcarc": load_pickle(gcarc_path)[gcmtid],
        "azimuth": load_pickle(azimuth_path)[gcmtid],
        "p": load_pickle(p_path)[gcmtid],
        "pp": load_pickle(pp_path)[gcmtid],
        "sp": load_pickle(sp_path)[gcmtid],
        "s": load_pickle(s_path)[gcmtid],
        "ss": load_pickle(ss_path)[gcmtid],
        "scs": load_pickle(scs_path)[gcmtid]
    }
    result = {}
    allkeys = list(load_pickle(gcarc_path)[gcmtid].keys())
    for each_key in allkeys:
        result[each_key] = {
            "gcarc": result_readin["gcarc"][each_key],
            "azimuth": result_readin["azimuth"][each_key],
            "p": result_readin["p"][each_key],
            "pp": result_readin["pp"][each_key],
            "sp": result_readin["sp"][each_key],
            "s": result_readin["s"][each_key],
            "ss": result_readin["ss"][each_key],
            "scs": result_readin["scs"][each_key]
        }
    evdp_dict = load_pickle(evdp_path)[gcmtid]
    rep_net_sta = list(evdp_dict.keys())[0]
    evdp = evdp_dict[rep_net_sta]
    return result, evdp


def build_to_plot_traces(obs_ds, syn_ds, trace_length, info_dir):
    # obs_ds,syn_ds opened asdf file
    # get keys
    key_obs = set(obs_ds.waveforms.list())
    key_syn = set(syn_ds.waveforms.list())
    keys = key_obs & key_syn
    result = {}
    # for each item in keys, get info
    # since the window is selected according to the two asdf files, we can just use keys
    gcmtid = obs_ds.events[0].resource_id.id.split("/")[-2]
    # ! fix a possible bug
    if (gcmtid == "smi:local"):
        gcmtid = obs_ds.events[0].resource_id.id.split("/")[-1].split("#")[0]
    info_dict, evdp = getinfo(info_dir, gcmtid)
    for key in keys:
        axkey = key.replace(".", "_")
        tag_obs = obs_ds.waveforms[key].get_waveform_tags()[0]
        tag_syn = syn_ds.waveforms[key].get_waveform_tags()[0]

        # here we use syn1_ds, which is not the normal case
        info = info_dict[key]
        obs_st = obs_ds.waveforms[key][tag_obs].copy()
        syn_st = syn_ds.waveforms[key][tag_syn].copy()

        # slice
        obs_st.trim(obs_st[0].stats.starttime,
                    obs_st[0].stats.starttime+trace_length)
        syn_st.trim(syn_st[0].stats.starttime,
                    syn_st[0].stats.starttime+trace_length)

        obs_r = obs_st.select(component="*R")[0]
        obs_t = obs_st.select(component="*T")[0]
        obs_z = obs_st.select(component="*Z")[0]
        syn_r = syn_st.select(component="*R")[0]
        syn_t = syn_st.select(component="*T")[0]
        syn_z = syn_st.select(component="*Z")[0]

        result[key] = to_plot_trace(
            obs_z, syn_z, obs_r, syn_r, obs_t, syn_t, info)
    return result, evdp


def build_plottting_structure(plot_traces, azimuth_width):
    # we assume 360%azimuth_width==0
    num_azimuths = 360//azimuth_width
    result = [[] for i in range(num_azimuths)]
    # for each item in plot_traces, seprate them into different []
    for key in plot_traces:
        value = plot_traces[key]
        info = value.info
        azimuth = info["azimuth"]
        index_azimuth = int(azimuth//azimuth_width)
        result[index_azimuth].append((key, value))

    # for each azimuth bin, sort them according to the gcarc
    def sort_func(item):
        return item[1].info["gcarc"]
    for index_azimuth in range(num_azimuths):
        result[index_azimuth] = sorted(result[index_azimuth], key=sort_func)
    return result


def kernel(obs_asdf, syn_asdf, azimuth_width, output_pdf, waves_perpage,
           trace_length, info_dir, misfit_windows_path, snr, cc, deltat, use_tqdm):
    obs_ds = pyasdf.ASDFDataSet(obs_asdf, mode="r", mpi=False)
    syn_ds = pyasdf.ASDFDataSet(syn_asdf, mode="r", mpi=False)
    windows_dict = load_pickle(misfit_windows_path)

    plot_traces, evdp = build_to_plot_traces(
        obs_ds, syn_ds, trace_length, info_dir)
    gcmtid = obs_ds.events[0].resource_id.id.split("/")[-2]
    plotting_structure = build_plottting_structure(plot_traces, azimuth_width)

    # estimate the total pages
    num_azimuths = 360 // azimuth_width
    count = 0
    if(use_tqdm):
        for index_azimuth in range(num_azimuths):
            azimuth_bin_plot_traces = plotting_structure[index_azimuth]
            num_azimuth_bin_plot_traces = len(azimuth_bin_plot_traces)
            # get num_pages for this azimuth bin
            if(num_azimuth_bin_plot_traces % waves_perpage == 0):
                num_pages = num_azimuth_bin_plot_traces // waves_perpage
            else:
                num_pages = (num_azimuth_bin_plot_traces // waves_perpage) + 1
            for ipage in range(num_pages):
                start_index = ipage*waves_perpage
                end_index = (ipage+1)*waves_perpage
                azimuth_bin_plot_traces_this_page = azimuth_bin_plot_traces[start_index:end_index]
                count += len(azimuth_bin_plot_traces_this_page)

    # plot figures
    pdf = matplotlib.backends.backend_pdf.PdfPages(output_pdf)
    if(use_tqdm):
        t = tqdm.tqdm(total=count)
    for index_azimuth in range(num_azimuths):
        # for each azimuth bin
        azimuth_bin_plot_traces = plotting_structure[index_azimuth]
        num_azimuth_bin_plot_traces = len(azimuth_bin_plot_traces)
        # get num_pages for this azimuth bin
        if(num_azimuth_bin_plot_traces % waves_perpage == 0):
            num_pages = num_azimuth_bin_plot_traces // waves_perpage
        else:
            num_pages = (num_azimuth_bin_plot_traces // waves_perpage)+1

        for ipage in range(num_pages):
            start_index = ipage*waves_perpage
            end_index = (ipage+1)*waves_perpage
            azimuth_bin_plot_traces_this_page = azimuth_bin_plot_traces[start_index:end_index]

            fig_z = plt.figure(figsize=(150, 150))
            fig_r = plt.figure(figsize=(150, 150))
            fig_t = plt.figure(figsize=(150, 150))
            index_count = 1
            axr, axz, axt = None, None, None  # get the last axes
            xticks = None
            for each_plot_trace_all in azimuth_bin_plot_traces_this_page:
                if(use_tqdm):
                    t.update()
                each_plot_trace = each_plot_trace_all[1]
                each_plot_id = each_plot_trace_all[0]
                # z
                axz = fig_z.add_subplot(waves_perpage, 1, index_count)
                obs = each_plot_trace.obs_z
                syn = each_plot_trace.syn_z
                x_obs = np.linspace(0, obs.stats.endtime -
                                    obs.stats.starttime, obs.stats.npts)
                x_syn = np.linspace(0, syn.stats.endtime -
                                    syn.stats.starttime, syn.stats.npts)
                y_obs = obs.data
                y_syn = syn.data
                axz.plot(x_obs, y_obs, color="k", linewidth=6.0)
                axz.plot(x_syn, y_syn, color="r", linewidth=6.0)

                axz.get_yaxis().set_ticklabels([])
                # r
                axr = fig_r.add_subplot(waves_perpage, 1, index_count)
                obs = each_plot_trace.obs_r
                syn = each_plot_trace.syn_r
                x_obs = np.linspace(0, obs.stats.endtime -
                                    obs.stats.starttime, obs.stats.npts)
                x_syn = np.linspace(0, syn.stats.endtime -
                                    syn.stats.starttime, syn.stats.npts)
                y_obs = obs.data
                y_syn = syn.data
                axr.plot(x_obs, y_obs, color="k", linewidth=6.0)
                axr.plot(x_syn, y_syn, color="r", linewidth=6.0)

                axr.get_yaxis().set_ticklabels([])
                # t
                axt = fig_t.add_subplot(waves_perpage, 1, index_count)
                obs = each_plot_trace.obs_t
                syn = each_plot_trace.syn_t
                x_obs = np.linspace(0, obs.stats.endtime -
                                    obs.stats.starttime, obs.stats.npts)
                x_syn = np.linspace(0, syn.stats.endtime -
                                    syn.stats.starttime, syn.stats.npts)
                y_obs = obs.data
                y_syn = syn.data
                axt.plot(x_obs, y_obs, color="k", linewidth=6.0)
                axt.plot(x_syn, y_syn, color="r", linewidth=6.0)

                axt.get_yaxis().set_ticklabels([])
                index_count += 1

                xticks = np.arange(np.min(x_obs), np.max(x_obs)+1, 100)
                axz.set_xticks(xticks)
                axr.set_xticks(xticks)
                axt.set_xticks(xticks)
                axz.tick_params(axis="x", labelsize=100)
                axr.tick_params(axis="x", labelsize=100)
                axt.tick_params(axis="x", labelsize=100)
                axz.tick_params(axis="y", labelsize=100)
                axr.tick_params(axis="y", labelsize=100)
                axt.tick_params(axis="y", labelsize=100)
                # plot texts
                axz.text(0.90, 0.53, f"id:{each_plot_id}\ngcarc:{each_plot_trace.info['gcarc']:.2f}\nazimuth:{each_plot_trace.info['azimuth']:.2f}",
                         transform=axz.transAxes, fontsize=100)
                axr.text(0.90, 0.53, f"id:{each_plot_id}\ngcarc:{each_plot_trace.info['gcarc']:.2f}\nazimuth:{each_plot_trace.info['azimuth']:.2f}",
                         transform=axr.transAxes, fontsize=100)
                axt.text(0.90, 0.53, f"id:{each_plot_id}\ngcarc:{each_plot_trace.info['gcarc']:.2f}\nazimuth:{each_plot_trace.info['azimuth']:.2f}",
                         transform=axt.transAxes, fontsize=100)

                # plot title
                if(index_count == 2):
                    axz.set_title(
                        f"azimuth:{azimuth_width*index_azimuth}-{azimuth_width*(index_azimuth+1)}\npage:{ipage}, vertical", fontsize=200)
                    axr.set_title(
                        f"azimuth:{azimuth_width*index_azimuth}-{azimuth_width*(index_azimuth+1)}\npage:{ipage}, radial", fontsize=200)
                    axt.set_title(
                        f"azimuth:{azimuth_width*index_azimuth}-{azimuth_width*(index_azimuth+1)}\npage:{ipage}, tagential", fontsize=200)

                # plot travel times
                info = each_plot_trace.info
                # z
                axz2 = axz.twiny()
                axz2.set_xlim(axz.get_xlim())
                twin_label = []
                plot_travel_times(axz2, "p", info["p"], np.max(
                    x_obs), "blue", twin_label)
                plot_travel_times(
                    axz2, "pp", info["pp"], np.max(x_obs), "y", twin_label)
                plot_travel_times(
                    axz2, "sp", info["sp"], np.max(x_obs), "r", twin_label)
                plot_travel_times(axz2, "s", info["s"], np.max(
                    x_obs), "green", twin_label)
                plot_travel_times(
                    axz2, "ss", info["ss"], np.max(x_obs), "black", twin_label)
                plot_windows(axz, windows_dict, each_plot_id,
                             "z", syn.stats.starttime, snr, cc, deltat, evdp)
                twin_x = [item[0] for item in twin_label]
                twin_name = [item[1] for item in twin_label]
                axz2.set_xticks(twin_x)
                axz2.set_xticklabels(twin_name, fontsize=50)
                axz2.xaxis.set_tick_params(rotation=90)
                # r
                axr2 = axr.twiny()
                axr2.set_xlim(axr.get_xlim())
                twin_label = []
                plot_travel_times(axr, "p", info["p"], np.max(
                    x_obs), "blue", twin_label)
                plot_travel_times(
                    axr, "pp", info["pp"], np.max(x_obs), "y", twin_label)
                plot_travel_times(
                    axr, "sp", info["sp"], np.max(x_obs), "r", twin_label)
                plot_travel_times(axr, "s", info["s"], np.max(
                    x_obs), "green", twin_label)
                plot_travel_times(
                    axr, "ss", info["ss"], np.max(x_obs), "black", twin_label)
                plot_windows(axr, windows_dict, each_plot_id,
                             "r", syn.stats.starttime, snr, cc, deltat, evdp)
                twin_x = [item[0] for item in twin_label]
                twin_name = [item[1] for item in twin_label]
                axr2.set_xticks(twin_x)
                axr2.set_xticklabels(twin_name, fontsize=50)
                axr2.xaxis.set_tick_params(rotation=90)
                # t
                axt2 = axt.twiny()
                axt2.set_xlim(axt.get_xlim())
                twin_label = []
                plot_travel_times(axt, "s", info["s"], np.max(
                    x_obs), "green", twin_label)
                plot_travel_times(
                    axt, "ss", info["ss"], np.max(x_obs), "black", twin_label)
                plot_travel_times(
                    axt, "scs", info["scs"], np.max(x_obs), "magenta", twin_label)
                plot_windows(axt, windows_dict, each_plot_id,
                             "t", syn.stats.starttime, snr, cc, deltat, evdp)
                twin_x = [item[0] for item in twin_label]
                twin_name = [item[1] for item in twin_label]
                axt2.set_xticks(twin_x)
                axt2.set_xticklabels(twin_name, fontsize=50)
                axt2.xaxis.set_tick_params(rotation=90)
            plt.subplots_adjust(wspace=0, hspace=0)
            fig_z.tight_layout()
            pdf.savefig(fig_z)
            fig_r.tight_layout()
            pdf.savefig(fig_r)
            fig_t.tight_layout()
            pdf.savefig(fig_t)
            plt.close(fig=fig_z)
            plt.close(fig=fig_r)
            plt.close(fig=fig_t)
    pdf.close()
    del obs_ds
    del syn_ds


@click.command()
@click.option('--obs_asdf', required=True, type=str)
@click.option('--syn_asdf', required=True, type=str)
@click.option('--azimuth_width', required=True, type=int)
@click.option('--output_pdf', required=True, type=str)
@click.option('--waves_perpage', required=True, type=int)
@click.option('--trace_length', required=True, type=int)
@click.option('--info_dir', required=True, type=str)
@click.option('--misfit_windows_path', required=True, type=str)
@click.option('--snr', required=True, type=float)
@click.option('--cc', required=True, type=float)
@click.option('--deltat', required=True, type=float)
@click.option('--use_tqdm', is_flag=True)
def main(obs_asdf, syn_asdf, azimuth_width, output_pdf, waves_perpage, trace_length, info_dir, misfit_windows_path, snr, cc, deltat, use_tqdm):
    kernel(obs_asdf, syn_asdf, azimuth_width, output_pdf, waves_perpage,
           trace_length, info_dir, misfit_windows_path, snr, cc, deltat, use_tqdm)


def plot_travel_times(ax, phasename, traveltime, length, thecolor, twin_label):
    if(traveltime == None):
        return
    if(traveltime < length):
        # ax.scatter(traveltime, 0, color=thecolor, label=phasename, s=9)
        ax.axvline(x=traveltime, color="g", linestyle='--', linewidth=3.0)
        # ax.text(traveltime, ylim[0]+(ylim[1]-ylim[0])*0.9, phasename)
        twin_label.append((traveltime, phasename))


def plot_windows(ax, windows_dict, net_sta, component, starttime, snr, cc, deltat, evdp):
    """
    plot the windows in the given ax.
    """
    # firstly we have to get the windows to plot.
    # ! note there is a possible bug that teh starttimes for the obs and syn are different due to the processing.
    windows_net_sta = windows_dict[net_sta]
    windows_to_plot_body = []
    windows_to_plot_surface = []
    if (component in windows_dict[net_sta]):
        for each_window in windows_dict[net_sta][component].windows:
            if((each_window.snr_energy >= snr) and (each_window.cc >= cc) and (np.abs(each_window.deltat) <= deltat)):
                windows_to_plot_body.append(each_window)
    if (("surface_" + component in windows_dict[net_sta]) and (evdp <= SURFACE_THRESHOLD)):
        for each_window in windows_dict[net_sta]["surface_" + component].windows:
            if((each_window.snr_energy >= snr) and (each_window.cc >= cc) and (np.abs(each_window.deltat) <= deltat)):
                windows_to_plot_surface.append(each_window)
    for each_window in windows_to_plot_body:
        left_index = each_window.left - starttime
        right_index = each_window.right - starttime
        ax.axvspan(left_index, right_index, alpha=0.10, color='green')
        # get plotting positions
        xmin, xmax = ax.get_xlim()
        x_text = (left_index - xmin) / (xmax - xmin)
        ax.text(x_text, 0.83, f"snr:{each_window.snr_energy:.2f}\ncc:{each_window.cc:.2f}\ndeltat:{each_window.deltat:.2f}\nsimilarity:{each_window.similarity:.2f}",
                transform=ax.transAxes, fontsize=25)
    for each_window in windows_to_plot_surface:
        left_index = each_window.left - starttime
        right_index = each_window.right - starttime
        ax.axvspan(left_index, right_index, alpha=0.10, color='yellow')
        # get plotting positions
        xmin, xmax = ax.get_xlim()
        x_text = (left_index - xmin) / (xmax - xmin)
        ax.text(x_text, 0.83, f"snr:{each_window.snr_energy:.2f}\ncc:{each_window.cc:.2f}\ndeltat:{each_window.deltat:.2f}\nsimilarity:{each_window.similarity:.2f}",
                transform=ax.transAxes, fontsize=25)


if __name__ == "__main__":
    main()
