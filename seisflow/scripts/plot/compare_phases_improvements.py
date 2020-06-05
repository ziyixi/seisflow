"""
compare_phases_improvements.py: for each phase, generate a pdf file to see the waveform improvements.
"""
from os.path import join, basename

import click
import matplotlib as mpl
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import numpy as np
import pyasdf
import tqdm

from ...setting import SURFACE_THRESHOLD
from ...utils.load_files import load_pickle
from .compare_two_asdf_with_windows import to_plot_trace, build_to_plot_traces
import obspy

label_size = 25
mpl.rcParams['xtick.labelsize'] = label_size

phases_zr = ["P", "pP", "sP"]
phases_zrt = ["S", "sS"]
phases_t = ["ScS"]

# map the phase name to the time name
name_mapper_time = {
    "P": "p",
    "pP": "pp",
    "sP": "sp",
    "S": "s",
    "sS": "ss",
    "ScS": "scs"
}
name_mapper_windows = {
    "P": ["P"],
    "pP": ["pP", "PP"],
    "sP": ["sP"],
    "S": ["S"],
    "sS": ["sS", "SS"],
    "ScS": ["ScS"]
}


def check_usable(windows_dict, phase_name, key, snr, cc, deltat, status_dict):
    """
    check if this net_sta in usable.
    """
    to_plot = False
    phase_name_list = name_mapper_windows[phase_name]
    windows_this_net_sta = windows_dict[key]
    for each_phase in phase_name_list:
        for each_component in windows_this_net_sta:
            thewindows = windows_this_net_sta[each_component].windows
            for each_window in thewindows:
                if (each_phase in each_window.phases):
                    if ((each_window.snr_energy >= snr) and (each_window.cc >= cc) and (np.abs(each_window.deltat) <= deltat)):
                        to_plot = True
                        status_dict[key][each_component] = True
    return to_plot


def slice_to_plot_traces(value_1, value_2, phase_name):
    """
    for all the components, we slice for the given phase_name
    """
    phase_name_used = name_mapper_time[phase_name]
    thetime = value_1.info[phase_name_used]
    if(thetime == None):
        return None, None
    eventtime = value_1.info["eventtime"]
    # "obs_z", "syn_z", "obs_r", "syn_r", "obs_t", "syn_t"
    sliced_value_1 = to_plot_trace(None, None, None, None, None, None, None)
    sliced_value_2 = to_plot_trace(None, None, None, None, None, None, None)
    sliced_value_1.info = value_1.info
    sliced_value_2.info = value_2.info
    sliced_value_1.obs_z = value_1.obs_z.slice(
        eventtime+thetime-20, eventtime+thetime+50)
    sliced_value_2.obs_z = value_2.obs_z.slice(
        eventtime+thetime-20, eventtime+thetime+50)
    sliced_value_1.syn_z = value_1.syn_z.slice(
        eventtime+thetime-20, eventtime+thetime+50)
    sliced_value_2.syn_z = value_2.syn_z.slice(
        eventtime + thetime - 20, eventtime + thetime + 50)
    # sliced_value_1.obs_z.taper(max_percentage=0.05)
    # sliced_value_1.obs_z.detrend()
    # sliced_value_2.obs_z.taper(max_percentage=0.05)
    # sliced_value_2.obs_z.detrend()
    normalize_st = obspy.Stream()+sliced_value_1.obs_z + \
        sliced_value_2.obs_z + sliced_value_1.syn_z + sliced_value_2.syn_z
    normalize_st.normalize(global_max=True)
    sliced_value_1.obs_r = value_1.obs_r.slice(
        eventtime+thetime-20, eventtime+thetime+50)
    sliced_value_1.syn_r = value_1.syn_r.slice(
        eventtime+thetime-20, eventtime+thetime+50)
    sliced_value_2.obs_r = value_2.obs_r.slice(
        eventtime+thetime-20, eventtime+thetime+50)
    sliced_value_2.syn_r = value_2.syn_r.slice(
        eventtime + thetime - 20, eventtime + thetime + 50)
    # sliced_value_1.obs_r.taper(max_percentage=0.05)
    # sliced_value_1.obs_r.detrend()
    # sliced_value_2.obs_r.taper(max_percentage=0.05)
    # sliced_value_2.obs_r.detrend()
    normalize_st = obspy.Stream()+sliced_value_1.obs_r + \
        sliced_value_2.obs_r + sliced_value_1.syn_r + sliced_value_2.syn_r
    normalize_st.normalize(global_max=True)
    sliced_value_1.obs_t = value_1.obs_t.slice(
        eventtime+thetime-20, eventtime+thetime+50)
    sliced_value_1.syn_t = value_1.syn_t.slice(
        eventtime+thetime-20, eventtime+thetime+50)
    sliced_value_2.obs_t = value_2.obs_t.slice(
        eventtime+thetime-20, eventtime+thetime+50)
    sliced_value_2.syn_t = value_2.syn_t.slice(
        eventtime + thetime - 20, eventtime + thetime + 50)
    # sliced_value_1.obs_t.taper(max_percentage=0.05)
    # sliced_value_1.obs_t.detrend()
    # sliced_value_2.obs_t.taper(max_percentage=0.05)
    # sliced_value_2.obs_t.detrend()
    normalize_st = obspy.Stream()+sliced_value_1.obs_t + \
        sliced_value_2.obs_t + sliced_value_1.syn_t + sliced_value_2.syn_t
    normalize_st.normalize(global_max=True)
    return sliced_value_1, sliced_value_2


def build_plottting_structure(plot_traces_1, plot_traces_2, azimuth_width, windows_dict, phase_name, snr, cc, deltat):
    # we assume 360%azimuth_width==0
    num_azimuths = 360//azimuth_width
    result = [[] for i in range(num_azimuths)]
    # for each item in plot_traces, seprate them into different []
    # since sync has the same net_sta list, we can use key directly from plot_traces_1
    # the status_dict used to see if this trace is usable
    status_dict = {}
    for key in plot_traces_1:
        status_dict[key] = {}
        if (check_usable(windows_dict, phase_name, key, snr, cc, deltat, status_dict) == False):
            continue
        value_1 = plot_traces_1[key]
        value_2 = plot_traces_2[key]
        # we should slice to the desired window
        sliced_value_1, sliced_value2 = slice_to_plot_traces(
            value_1, value_2, phase_name)
        if((sliced_value_1 == None) and (sliced_value2 == None)):
            continue
        info = value_1.info
        azimuth = info["azimuth"]
        index_azimuth = int(azimuth//azimuth_width)
        result[index_azimuth].append((key, sliced_value_1, sliced_value2))

    # for each azimuth bin, sort them according to the gcarc
    def sort_func(item):
        return item[1].info["gcarc"]
    for index_azimuth in range(num_azimuths):
        result[index_azimuth] = sorted(result[index_azimuth], key=sort_func)
    return result, status_dict


def kernel(obs_asdf, syn_asdf_1, syn_asdf_2, azimuth_width, output_dir, waves_perpage,
           trace_length, info_dir, misfit_windows_path_1, misfit_windows_path_2, snr, cc, deltat, use_tqdm):
    """
    for each phase, find if the phase window is avaliable or not, and plot them into a pdf file.
    """
    obs_ds = pyasdf.ASDFDataSet(obs_asdf, mode="r", mpi=False)
    syn_ds_1 = pyasdf.ASDFDataSet(syn_asdf_1, mode="r", mpi=False)
    syn_ds_2 = pyasdf.ASDFDataSet(syn_asdf_2, mode="r", mpi=False)
    windows_dict_1 = load_pickle(misfit_windows_path_1)
    windows_dict_2 = load_pickle(misfit_windows_path_2)
    # get pdf base name
    pdf_base = ".".join(basename(obs_asdf).split(".")[:2])

    plot_traces_1, evdp = build_to_plot_traces(
        obs_ds, syn_ds_1, trace_length, info_dir)
    plot_traces_2, evdp = build_to_plot_traces(
        obs_ds, syn_ds_2, trace_length, info_dir)

    num_azimuths = 360 // azimuth_width

    # * firstly we generate all the plotting_structure and status_dict
    plotting_structure_dict = {}
    status_dict_dict = {}
    for each_phase in (phases_zr + phases_zrt + phases_t):
        plotting_structure, status_dict = build_plottting_structure(
            plot_traces_1, plot_traces_2, azimuth_width, windows_dict_1, each_phase, snr, cc, deltat)
        plotting_structure_dict[each_phase] = plotting_structure
        status_dict_dict[each_phase] = status_dict

    # * get the tqdm count
    count = 0
    if (use_tqdm):
        for each_phase in (phases_zr + phases_zrt + phases_t):
            plotting_structure, status_dict = plotting_structure_dict[
                each_phase], status_dict_dict[each_phase]
            for index_azimuth in range(num_azimuths):
                azimuth_bin_plot_traces = plotting_structure[index_azimuth]
                num_azimuth_bin_plot_traces = len(azimuth_bin_plot_traces)
                # get num_pages for this azimuth bin
                if(num_azimuth_bin_plot_traces % waves_perpage == 0):
                    num_pages = num_azimuth_bin_plot_traces // waves_perpage
                else:
                    num_pages = (num_azimuth_bin_plot_traces //
                                 waves_perpage)+1
                for ipage in range(num_pages):
                    start_index = ipage*waves_perpage
                    end_index = (ipage+1)*waves_perpage
                    azimuth_bin_plot_traces_this_page = azimuth_bin_plot_traces[start_index:end_index]
                    for index, each_plot_trace_all in enumerate(azimuth_bin_plot_traces_this_page):
                        count += 1
    if(use_tqdm):
        t = tqdm.tqdm(total=count)

    # * for each phase in phase_zr, we plot the comparision.
    for each_phase in (phases_zr+phases_zrt+phases_t):
        plotting_structure, status_dict = plotting_structure_dict[
            each_phase], status_dict_dict[each_phase]
        output_path = join(output_dir, f"{pdf_base}.{each_phase}.pdf")
        pdf = matplotlib.backends.backend_pdf.PdfPages(output_path)
        for index_azimuth in range(num_azimuths):
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

                fig = plt.figure(figsize=(150, 150))
                xticks = None
                for index, each_plot_trace_all in enumerate(azimuth_bin_plot_traces_this_page):
                    if(use_tqdm):
                        t.update()
                    if (index == 0):
                        ax_title = fig.add_subplot(
                            waves_perpage, 1, 1)
                        ax_title.set_title(
                            f"azimuth:{azimuth_width*index_azimuth}-{azimuth_width*(index_azimuth+1)}\npage:{ipage}", fontsize=200)
                    # * we should plot 4 subplots per row (6 figures spacing)
                    each_plot_id = each_plot_trace_all[0]
                    each_plot_trace_1 = each_plot_trace_all[1]
                    each_plot_trace_2 = each_plot_trace_all[2]

                    # * firstly we plot the z component
                    ax_1 = fig.add_subplot(waves_perpage, 6, index * 6 + 1)
                    ax_2 = fig.add_subplot(
                        waves_perpage, 6, index * 6 + 2, sharey=ax_1)
                    # * get the traces
                    obs = each_plot_trace_1.obs_z
                    syn_1 = each_plot_trace_1.syn_z
                    syn_2 = each_plot_trace_2.syn_z
                    # * now we should get x
                    x_obs = np.linspace(-20, obs.stats.endtime -
                                        obs.stats.starttime-20, obs.stats.npts)
                    x_syn_1 = np.linspace(-20, syn_1.stats.endtime -
                                          syn_1.stats.starttime-20, syn_1.stats.npts)
                    x_syn_2 = np.linspace(-20, syn_2.stats.endtime -
                                          syn_2.stats.starttime-20, syn_2.stats.npts)
                    y_obs = obs.data
                    y_syn_1 = syn_1.data
                    y_syn_2 = syn_2.data
                    ax_1.plot(x_obs, y_obs, color="k", linewidth=6.0)
                    ax_1.plot(x_syn_1, y_syn_1, color="b", linewidth=6.0)
                    ax_2.plot(x_obs, y_obs, color="k", linewidth=6.0)
                    ax_2.plot(x_syn_2, y_syn_2, color="r", linewidth=6.0)
                    ax_1.axvline(x=0, color="g",
                                 linestyle='--', linewidth=3.0)
                    ax_2.axvline(x=0, color="g",
                                 linestyle='--', linewidth=3.0)
                    windows_list_1 = windows_dict_1[each_plot_id]["z"].windows
                    windows_list_2 = windows_dict_2[each_plot_id]["z"].windows
                    each_window_1 = None
                    for item in windows_list_1:
                        if (each_phase in item.phases):
                            each_window_1 = item
                            break
                    each_window_2 = None
                    for item in windows_list_2:
                        if (each_phase in item.phases):
                            each_window_2 = item
                            break
                    xmin, xmax = ax_1.get_xlim()
                    x_text = (-20 - xmin) / (xmax - xmin)
                    if(each_window_1 != None):
                        ax_1.text(x_text, 0.73, f"snr:{each_window_1.snr_energy:.2f}\ncc:{each_window_1.cc:.2f}\ndeltat:{each_window_1.deltat:.2f}\nsimilarity:{each_window_1.similarity:.2f}",
                                  transform=ax_1.transAxes, fontsize=50)
                    else:
                        ax_1.text(x_text, 0.73, f"snr:None\ncc:None\ndeltat:None\nsimilarity:None",
                                  transform=ax_1.transAxes, fontsize=50)
                    ax_1.text((35 - xmin) / (xmax - xmin), 0.73, f"component:Z\nid:{each_plot_id}\ngcarc:{each_plot_trace_1.info['gcarc']:.2f}\nazimuth:{each_plot_trace_1.info['azimuth']:.2f}",
                              transform=ax_1.transAxes, fontsize=50)
                    xmin, xmax = ax_2.get_xlim()
                    x_text = (-20 - xmin) / (xmax - xmin)
                    if(each_window_2 != None):
                        ax_2.text(x_text, 0.73, f"snr:{each_window_2.snr_energy:.2f}\ncc:{each_window_2.cc:.2f}\ndeltat:{each_window_2.deltat:.2f}\nsimilarity:{each_window_2.similarity:.2f}",
                                  transform=ax_2.transAxes, fontsize=50)
                    else:
                        ax_2.text(x_text, 0.73, f"snr:None\ncc:None\ndeltat:None\nsimilarity:None",
                                  transform=ax_2.transAxes, fontsize=50)
                    ax_2.text((35 - xmin) / (xmax - xmin), 0.73, f"component:Z\nid:{each_plot_id}\ngcarc:{each_plot_trace_2.info['gcarc']:.2f}\nazimuth:{each_plot_trace_2.info['azimuth']:.2f}",
                              transform=ax_2.transAxes, fontsize=50)
                    if ("z" not in status_dict[each_plot_id]):
                        ax_1.axvspan(xmin, xmax,
                                     alpha=0.10, color='red')
                        ax_2.axvspan(xmin, xmax,
                                     alpha=0.10, color='red')
                    else:
                        ax_1.axvspan(xmin, xmax,
                                     alpha=0.10, color='green')
                        ax_2.axvspan(xmin, xmax,
                                     alpha=0.10, color='green')
                    ax_1.get_yaxis().set_ticklabels([])
                    ax_2.get_yaxis().set_ticklabels([])
                    ax_1.spines["top"].set_visible(False)
                    ax_1.spines["right"].set_visible(False)
                    ax_1.spines["bottom"].set_visible(False)
                    ax_1.spines["left"].set_visible(False)
                    ax_2.spines["top"].set_visible(False)
                    ax_2.spines["right"].set_visible(False)
                    ax_2.spines["bottom"].set_visible(False)
                    ax_2.spines["left"].set_visible(False)
                    xticks = np.arange(np.min(x_obs), np.max(x_obs) + 1, 10)
                    ax_1.set_xticks(xticks)
                    ax_2.set_xticks(xticks)
                    ax_1.tick_params(axis="x", labelsize=100)
                    ax_2.tick_params(axis="x", labelsize=100)

                    # * and we plot the r component
                    ax_1 = fig.add_subplot(waves_perpage, 6, index * 6 + 3)
                    ax_2 = fig.add_subplot(
                        waves_perpage, 6, index * 6 + 4, sharey=ax_1)
                    # * get the traces
                    obs = each_plot_trace_1.obs_r
                    syn_1 = each_plot_trace_1.syn_r
                    syn_2 = each_plot_trace_2.syn_r
                    # * now we should get x
                    x_obs = np.linspace(-20, obs.stats.endtime -
                                        obs.stats.starttime-20, obs.stats.npts)
                    x_syn_1 = np.linspace(-20, syn_1.stats.endtime -
                                          syn_1.stats.starttime-20, syn_1.stats.npts)
                    x_syn_2 = np.linspace(-20, syn_2.stats.endtime -
                                          syn_2.stats.starttime-20, syn_2.stats.npts)
                    y_obs = obs.data
                    y_syn_1 = syn_1.data
                    y_syn_2 = syn_2.data
                    ax_1.plot(x_obs, y_obs, color="k", linewidth=6.0)
                    ax_1.plot(x_syn_1, y_syn_1, color="b", linewidth=6.0)
                    ax_2.plot(x_obs, y_obs, color="k", linewidth=6.0)
                    ax_2.plot(x_syn_2, y_syn_2, color="r", linewidth=6.0)
                    ax_1.axvline(x=0, color="g",
                                 linestyle='--', linewidth=3.0)
                    ax_2.axvline(x=0, color="g",
                                 linestyle='--', linewidth=3.0)
                    windows_list_1 = windows_dict_1[each_plot_id]["r"].windows
                    windows_list_2 = windows_dict_2[each_plot_id]["r"].windows
                    each_window_1 = None
                    for item in windows_list_1:
                        if (each_phase in item.phases):
                            each_window_1 = item
                            break
                    each_window_2 = None
                    for item in windows_list_2:
                        if (each_phase in item.phases):
                            each_window_2 = item
                            break
                    xmin, xmax = ax_1.get_xlim()
                    x_text = (-20 - xmin) / (xmax - xmin)
                    if(each_window_1 != None):
                        ax_1.text(x_text, 0.73, f"snr:{each_window_1.snr_energy:.2f}\ncc:{each_window_1.cc:.2f}\ndeltat:{each_window_1.deltat:.2f}\nsimilarity:{each_window_1.similarity:.2f}",
                                  transform=ax_1.transAxes, fontsize=50)
                    else:
                        ax_1.text(x_text, 0.73, f"snr:None\ncc:None\ndeltat:None\nsimilarity:None",
                                  transform=ax_1.transAxes, fontsize=50)
                    ax_1.text((35 - xmin) / (xmax - xmin), 0.73, f"component:R\nid:{each_plot_id}\ngcarc:{each_plot_trace_1.info['gcarc']:.2f}\nazimuth:{each_plot_trace_1.info['azimuth']:.2f}",
                              transform=ax_1.transAxes, fontsize=50)
                    xmin, xmax = ax_2.get_xlim()
                    x_text = (-20 - xmin) / (xmax - xmin)
                    if(each_window_2 != None):
                        ax_2.text(x_text, 0.73, f"snr:{each_window_2.snr_energy:.2f}\ncc:{each_window_2.cc:.2f}\ndeltat:{each_window_2.deltat:.2f}\nsimilarity:{each_window_2.similarity:.2f}",
                                  transform=ax_2.transAxes, fontsize=50)
                    else:
                        ax_2.text(x_text, 0.73, f"snr:None\ncc:None\ndeltat:None\nsimilarity:None",
                                  transform=ax_2.transAxes, fontsize=50)
                    ax_2.text((35 - xmin) / (xmax - xmin), 0.73, f"component:R\nid:{each_plot_id}\ngcarc:{each_plot_trace_2.info['gcarc']:.2f}\nazimuth:{each_plot_trace_2.info['azimuth']:.2f}",
                              transform=ax_2.transAxes, fontsize=50)
                    if ("r" not in status_dict[each_plot_id]):
                        ax_1.axvspan(xmin, xmax,
                                     alpha=0.10, color='red')
                        ax_2.axvspan(xmin, xmax,
                                     alpha=0.10, color='red')
                    else:
                        ax_1.axvspan(xmin, xmax,
                                     alpha=0.10, color='green')
                        ax_2.axvspan(xmin, xmax,
                                     alpha=0.10, color='green')
                    ax_1.get_yaxis().set_ticklabels([])
                    ax_2.get_yaxis().set_ticklabels([])
                    ax_1.spines["top"].set_visible(False)
                    ax_1.spines["right"].set_visible(False)
                    ax_1.spines["bottom"].set_visible(False)
                    ax_1.spines["left"].set_visible(False)
                    ax_2.spines["top"].set_visible(False)
                    ax_2.spines["right"].set_visible(False)
                    ax_2.spines["bottom"].set_visible(False)
                    ax_2.spines["left"].set_visible(False)
                    xticks = np.arange(np.min(x_obs), np.max(x_obs) + 1, 10)
                    ax_1.set_xticks(xticks)
                    ax_2.set_xticks(xticks)
                    ax_1.tick_params(axis="x", labelsize=100)
                    ax_2.tick_params(axis="x", labelsize=100)

                    # * and we plot the t component
                    ax_1 = fig.add_subplot(waves_perpage, 6, index * 6 + 5)
                    ax_2 = fig.add_subplot(
                        waves_perpage, 6, index * 6 + 6, sharey=ax_1)
                    # * get the traces
                    obs = each_plot_trace_1.obs_t
                    syn_1 = each_plot_trace_1.syn_t
                    syn_2 = each_plot_trace_2.syn_t
                    # * now we should get x
                    x_obs = np.linspace(-20, obs.stats.endtime -
                                        obs.stats.starttime-20, obs.stats.npts)
                    x_syn_1 = np.linspace(-20, syn_1.stats.endtime -
                                          syn_1.stats.starttime-20, syn_1.stats.npts)
                    x_syn_2 = np.linspace(-20, syn_2.stats.endtime -
                                          syn_2.stats.starttime-20, syn_2.stats.npts)
                    y_obs = obs.data
                    y_syn_1 = syn_1.data
                    y_syn_2 = syn_2.data
                    ax_1.plot(x_obs, y_obs, color="k", linewidth=6.0)
                    ax_1.plot(x_syn_1, y_syn_1, color="b", linewidth=6.0)
                    ax_2.plot(x_obs, y_obs, color="k", linewidth=6.0)
                    ax_2.plot(x_syn_2, y_syn_2, color="r", linewidth=6.0)
                    ax_1.axvline(x=0, color="g",
                                 linestyle='--', linewidth=3.0)
                    ax_2.axvline(x=0, color="g",
                                 linestyle='--', linewidth=3.0)
                    windows_list_1 = windows_dict_1[each_plot_id]["t"].windows
                    windows_list_2 = windows_dict_2[each_plot_id]["t"].windows
                    each_window_1 = None
                    for item in windows_list_1:
                        if (each_phase in item.phases):
                            each_window_1 = item
                            break
                    each_window_2 = None
                    for item in windows_list_2:
                        if (each_phase in item.phases):
                            each_window_2 = item
                            break
                    xmin, xmax = ax_1.get_xlim()
                    x_text = (-20 - xmin) / (xmax - xmin)
                    if(each_window_1 != None):
                        ax_1.text(x_text, 0.73, f"snr:{each_window_1.snr_energy:.2f}\ncc:{each_window_1.cc:.2f}\ndeltat:{each_window_1.deltat:.2f}\nsimilarity:{each_window_1.similarity:.2f}",
                                  transform=ax_1.transAxes, fontsize=50)
                    else:
                        ax_1.text(x_text, 0.73, f"snr:None\ncc:None\ndeltat:None\nsimilarity:None",
                                  transform=ax_1.transAxes, fontsize=50)
                    ax_1.text((35 - xmin) / (xmax - xmin), 0.73, f"component:T\nid:{each_plot_id}\ngcarc:{each_plot_trace_1.info['gcarc']:.2f}\nazimuth:{each_plot_trace_1.info['azimuth']:.2f}",
                              transform=ax_1.transAxes, fontsize=50)
                    xmin, xmax = ax_2.get_xlim()
                    x_text = (-20 - xmin) / (xmax - xmin)
                    if(each_window_2 != None):
                        ax_2.text(x_text, 0.73, f"snr:{each_window_2.snr_energy:.2f}\ncc:{each_window_2.cc:.2f}\ndeltat:{each_window_2.deltat:.2f}\nsimilarity:{each_window_2.similarity:.2f}",
                                  transform=ax_2.transAxes, fontsize=50)
                    else:
                        ax_2.text(x_text, 0.73, f"snr:None\ncc:None\ndeltat:None\nsimilarity:None",
                                  transform=ax_2.transAxes, fontsize=50)
                    ax_2.text((35 - xmin) / (xmax - xmin), 0.73, f"component:T\nid:{each_plot_id}\ngcarc:{each_plot_trace_2.info['gcarc']:.2f}\nazimuth:{each_plot_trace_2.info['azimuth']:.2f}",
                              transform=ax_2.transAxes, fontsize=50)
                    if ("t" not in status_dict[each_plot_id]):
                        ax_1.axvspan(xmin, xmax,
                                     alpha=0.10, color='red')
                        ax_2.axvspan(xmin, xmax,
                                     alpha=0.10, color='red')
                    else:
                        ax_1.axvspan(xmin, xmax,
                                     alpha=0.10, color='green')
                        ax_2.axvspan(xmin, xmax,
                                     alpha=0.10, color='green')
                    ax_1.get_yaxis().set_ticklabels([])
                    ax_2.get_yaxis().set_ticklabels([])
                    ax_1.spines["top"].set_visible(False)
                    ax_1.spines["right"].set_visible(False)
                    ax_1.spines["bottom"].set_visible(False)
                    ax_1.spines["left"].set_visible(False)
                    ax_2.spines["top"].set_visible(False)
                    ax_2.spines["right"].set_visible(False)
                    ax_2.spines["bottom"].set_visible(False)
                    ax_2.spines["left"].set_visible(False)
                    xticks = np.arange(np.min(x_obs), np.max(x_obs) + 1, 10)
                    ax_1.set_xticks(xticks)
                    ax_2.set_xticks(xticks)
                    ax_1.tick_params(axis="x", labelsize=100)
                    ax_2.tick_params(axis="x", labelsize=100)

                plt.subplots_adjust(wspace=0, hspace=0, top=0.85)
                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig=fig)
        pdf.close()
    del obs_ds
    del syn_ds_1
    del syn_ds_2
    if (use_tqdm):
        t.close()


if __name__ == "__main__":
    @click.command()
    @click.option('--obs_asdf', required=True, type=str)
    @click.option('--syn_asdf_1', required=True, type=str)
    @click.option('--syn_asdf_2', required=True, type=str)
    @click.option('--azimuth_width', required=True, type=int)
    @click.option('--output_dir', required=True, type=str)
    @click.option('--waves_perpage', required=True, type=int)
    @click.option('--trace_length', required=True, type=int)
    @click.option('--info_dir', required=True, type=str)
    @click.option('--misfit_windows_path_1', required=True, type=str)
    @click.option('--misfit_windows_path_2', required=True, type=str)
    @click.option('--snr', required=True, type=float)
    @click.option('--cc', required=True, type=float)
    @click.option('--deltat', required=True, type=float)
    @click.option('--use_tqdm', is_flag=True)
    def main(obs_asdf, syn_asdf_1, syn_asdf_2, azimuth_width, output_dir, waves_perpage,
             trace_length, info_dir, misfit_windows_path_1, misfit_windows_path_2, snr, cc, deltat, use_tqdm):
        kernel(obs_asdf, syn_asdf_1, syn_asdf_2, azimuth_width, output_dir, waves_perpage,
               trace_length, info_dir, misfit_windows_path_1, misfit_windows_path_2, snr, cc, deltat, use_tqdm)

    main()  # pylint: disable=no-value-for-parameter
