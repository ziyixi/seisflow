"""
plot_source_inversion_improvements.py: plot the source inversion improvements for different kinds of phases.
"""
from glob import glob
from os.path import basename, join

import matplotlib.pyplot as plt
import numpy as np
import obspy
import seaborn as sns
import tqdm
from matplotlib.backends.backend_pdf import PdfPages
from obspy.imaging.beachball import beach

from ...utils.data_analysis import (get_windows_cc, get_windows_deltat,
                                    get_windows_net_sta,
                                    get_windows_similarity,
                                    get_windows_snr_energy)
from ...utils.load_files import load_first_arrival_baz_evdp, load_pickle
from ...utils.setting import SURFACE_THRESHOLD


def load_misfit_windows(misfit_windows_base, iterations_list):
    """
    load all the misfit windows in the misfit_windows_base given the iteration numbers.
    """
    misfit_windows_collection = {}
    for each_iteration in iterations_list:
        iteration_directory = join(
            misfit_windows_base, f"iter{each_iteration}_misfit_windows")
        misfit_windows_collection[each_iteration] = {}
        all_files = sorted(glob(join(iteration_directory, "*pkl")))
        for each_file in all_files:
            gcmtid = basename(each_file).split(".")[0]
            misfit_windows_collection[each_iteration][gcmtid] = load_pickle(
                each_file)
    return misfit_windows_collection


def get_plotting_data(each_misfit_windows_collection, iterations_list, snr_threshold, event_depth_dict):
    """
    get the list to plot in seaborn. each page is for a category. In each page, row is deltat,similarity,cc; column is for different phases
    """
    result = {}
    phases_zr = ["P", "pP", "sP", "PP", "S", "sS", "SS"]
    phases_t = ["ScS", "S", "sS", "SS"]
    conditions = {
        "P": {
            "exclude_p": False,
            "exclude_s": True
        },
        "pP": {
            "exclude_p": True,
            "exclude_s": True
        },
        "sP": {
            "exclude_p": True,
            "exclude_s": True
        },
        "PP": {
            "exclude_p": True,
            "exclude_s": True
        },
        "S": {
            "exclude_p": False,
            "exclude_s": False
        },
        "sS": {
            "exclude_p": True,
            "exclude_s": True
        },
        "SS": {
            "exclude_p": True,
            "exclude_s": True
        },
        "ScS": {
            "exclude_p": True,
            "exclude_s": True
        },
        "surface_z": {
            "exclude_p": False,
            "exclude_s": False
        },
        "surface_r": {
            "exclude_p": False,
            "exclude_s": False
        },
        "surface_t": {
            "exclude_p": False,
            "exclude_s": False
        },
    }
    # we can exrtact the information from the misfit_windows in the order of the pdf output.
    # order will be z,r,t[,surface_z,surface_r,surface_t]
    rep_net_sta = sorted(event_depth_dict.keys())[0]
    event_depth_this_event = event_depth_dict[rep_net_sta]
    if (event_depth_this_event > SURFACE_THRESHOLD):
        category_list = ["z", "r", "t"]
        category_phases = [phases_zr, phases_zr, phases_t]
    else:
        category_list = ["z", "r", "t", "surface_z", "surface_r", "surface_t"]
        category_phases = [phases_zr, phases_zr, phases_t,
                           ["surface_z"], ["surface_r"], ["surface_t"]]
    for each_iteration in iterations_list:
        result[each_iteration] = {}
        for each_category, each_category_phases in zip(category_list, category_phases):
            result[each_iteration][each_category] = []
            for each_category_phase in each_category_phases:
                phase_condition = conditions[each_category_phase]
                cc = get_windows_cc(
                    each_misfit_windows_collection[each_iteration], phase_condition[
                        "exclude_p"], phase_condition["exclude_s"],
                    each_category, snr_threshold, each_category_phase)
                cc = cc[cc >= 0]
                deltat = get_windows_deltat(
                    each_misfit_windows_collection[each_iteration], phase_condition[
                        "exclude_p"], phase_condition["exclude_s"],
                    each_category, snr_threshold, each_category_phase)
                deltat = deltat[np.abs(deltat) <= 10]
                similarity = get_windows_similarity(
                    each_misfit_windows_collection[each_iteration], phase_condition[
                        "exclude_p"], phase_condition["exclude_s"],
                    each_category, snr_threshold, each_category_phase)
                similarity = similarity[similarity >= 0]
                result[each_iteration][each_category].append(
                    {"net_sta": get_windows_net_sta(
                        each_misfit_windows_collection[each_iteration], phase_condition[
                            "exclude_p"], phase_condition["exclude_s"],
                        each_category, snr_threshold, each_category_phase),
                     "cc": cc,
                     "deltat": deltat,
                     "similarity": similarity,
                     }
                )
        # result:dict->each_iteration:dict->each_category:list as the dict showed before, we should return the category_phases
        # we should combine the surface wave phases to one page
    if (len(category_phases) == 6):
        for each_iteration in iterations_list:
            category_phases = [phases_zr, phases_zr, phases_t,
                               ["surface_z", "surface_r", "surface_t"]]
            category_list = ["z", "r", "t", "surface"]
            result[each_iteration]["surface"] = []
            result[each_iteration]["surface"].append(
                result[each_iteration]["surface_z"][0])
            result[each_iteration]["surface"].append(
                result[each_iteration]["surface_r"][0])
            result[each_iteration]["surface"].append(
                result[each_iteration]["surface_t"][0])
            del result[each_iteration]["surface_z"]
            del result[each_iteration]["surface_r"]
            del result[each_iteration]["surface_t"]

    return result, category_phases, category_list


def plot_source_parameters(event_gcmtid, pdf, cmts_directory, iterations_list):
    """
    plot the beachball and the source parameters.
    """
    figs = plt.figure(figsize=(50, 50))
    # * firstly we prepare the information from cmt solution
    # our desired columns will be for each iteration
    # our desired rows will be [time, half duration, latitude, longitude, depth.]
    # and we also plot the beachball in the bottom.
    table_data = np.zeros((5, len(iterations_list)), dtype=np.object)
    for column_index, each_iteration in enumerate(iterations_list):
        if (int(each_iteration) == 1):
            cmt_fpath = join(cmts_directory, "raw_cmtsolutions", event_gcmtid)
        else:
            cmt_fpath = join(
                cmts_directory, f"iter{int(each_iteration)-1}_next_cmtsolutions", event_gcmtid)
        cmtsolution = obspy.read_events(cmt_fpath)[0]
        table_data[0, column_index] = str(cmtsolution.preferred_origin().time)
        table_data[1, column_index] = cmtsolution.focal_mechanisms[0].moment_tensor.source_time_function.duration / 2
        table_data[2, column_index] = cmtsolution.preferred_origin().latitude
        table_data[3, column_index] = cmtsolution.preferred_origin().longitude
        table_data[4, column_index] = cmtsolution.preferred_origin().depth / \
            1000.0

    ax = figs.add_subplot(2, 1, 1)
    figs.patch.set_visible(False)
    ax.axis('off')
    thetable = ax.table(cellText=table_data, rowLabels=["time", "half duration(s)", "latitude(°)", "longitude(°)", "depth(km)"],
                        colLabels=[
                            f"{event_gcmtid}\niteration {item}" for item in iterations_list],
                        cellLoc='center', rowLoc='center', loc="center")
    thetable.set_fontsize(75)
    table_props = thetable.properties()
    table_cells = table_props['child_artists']
    for cell in table_cells:
        cell.set_height(0.3)
    for column_index, each_iteration in enumerate(iterations_list):
        if (int(each_iteration) == 1):
            cmt_fpath = join(cmts_directory, "raw_cmtsolutions", event_gcmtid)
        else:
            cmt_fpath = join(
                cmts_directory, f"iter{int(each_iteration)-1}_next_cmtsolutions", event_gcmtid)
        cmtsolution = obspy.read_events(cmt_fpath)[0]
        # plot the beachball
        ax = figs.add_subplot(6, len(iterations_list),
                              column_index + 4 * len(iterations_list) + 1)
        ax.axis('off')
        tensor = cmtsolution.focal_mechanisms[0].moment_tensor.tensor
        to_plot_beachball = beach(
            [tensor.m_rr, tensor.m_tt, tensor.m_pp, tensor.m_rt, tensor.m_rp, tensor.m_tp], width=.5, xy=(0.5, 0.5))
        ax.add_collection(to_plot_beachball)
        ax.set_xlim((0, 1))
        ax.set_ylim((0, 1))
        ax.axis('equal')
    figs.tight_layout()
    pdf.savefig(figs, bbox_inches="tight")
    plt.close(fig=figs)


def plot_to_pdf(pdf_fname, cmts_directory, misfit_windows_collection, iterations_list, snr_threshold, event_depth):
    """
    figs = plt.figure(figsize=(50, 50))
    """
    rep_key = sorted(misfit_windows_collection.keys())[0]
    all_events = sorted(misfit_windows_collection[rep_key].keys())
    with PdfPages(pdf_fname) as pdf:
        for each_event in tqdm.tqdm(all_events):
            # we should plot the beachball and plot the source parameter table here
            plot_source_parameters(
                each_event, pdf, cmts_directory, iterations_list)
            # prepare information to plot
            each_misfit_windows_collection = {}
            for each_iteration in iterations_list:
                each_misfit_windows_collection[each_iteration] = (
                    misfit_windows_collection[each_iteration][each_event])
            event_depth_dict = event_depth[each_event]
            data_collection, category_phases, category_list = get_plotting_data(
                each_misfit_windows_collection, iterations_list, snr_threshold, event_depth_dict)
            for each_category, phase_list_for_each_category in zip(category_list, category_phases):
                # one page for each category
                figs = plt.figure(figsize=(50, 50))
                collecction_all = {}
                if (each_category != "surface"):
                    collecction_all["deltat"] = [np.array([], dtype=np.float)
                                                 for i in range(len(iterations_list))]
                    collecction_all["similarity"] = [np.array([], dtype=np.float)
                                                     for i in range(len(iterations_list))]
                    collecction_all["cc"] = [np.array([], dtype=np.float)
                                             for i in range(len(iterations_list))]
                # we plot for each phases
                for row_index, each_phase in enumerate(phase_list_for_each_category):
                    # we plot for deltat,similarity,cc
                    for column_index, plot_type in enumerate(["deltat", "similarity", "cc"]):
                        # num must be 1 <= num <= num_max, not 0
                        # keep different category's figsize the same
                        ax = figs.add_subplot(
                            8, 3, row_index * 3 + column_index+1)

                        for interation_index, each_iteration in enumerate(iterations_list):
                            sns.distplot(data_collection[each_iteration][each_category][row_index]
                                         [plot_type], ax=ax, hist=False, label=f"before iteration {each_iteration}",
                                         kde_kws={"linewidth": 6})
                            # collect to the category summary
                            if(each_category != "surface"):
                                if (column_index == 0):
                                    collecction_all["deltat"][interation_index] = np.concatenate(
                                        (collecction_all["deltat"][interation_index], data_collection[each_iteration][each_category][row_index]
                                         [plot_type]))
                                elif (column_index == 1):
                                    collecction_all["similarity"][interation_index] = np.concatenate(
                                        (collecction_all["similarity"][interation_index], data_collection[each_iteration][each_category][row_index]
                                         [plot_type]))
                                elif (column_index == 2):
                                    collecction_all["cc"][interation_index] = np.concatenate(
                                        (collecction_all["cc"][interation_index], data_collection[each_iteration][each_category][row_index]
                                         [plot_type]))
                        if (plot_type == "deltat"):
                            ax.set_xlim((-10, 10))
                        elif(plot_type == "similarity"):
                            ax.set_xlim((0, 1))
                        elif(plot_type == "cc"):
                            ax.set_xlim((0, 1))
                        # ax.legend()
                        if (column_index == 0):
                            ax.get_yaxis().set_ticklabels([])
                            ax.set_ylabel(each_phase, fontsize=50, rotation=90)
                        else:
                            ax.get_yaxis().set_ticklabels([])
                        ax.tick_params(axis="x", labelsize=30)
                        if(plot_type != "similarity"):
                            ax.set_xlabel(plot_type, fontsize=30)
                        else:
                            ax.set_xlabel("zero-lag cc", fontsize=30)
                        if (row_index == 0 and column_index == 1):
                            ax.set_title(
                                f"gcmtid: {each_event}\ncategory: {each_category}", fontsize=50)
                if (each_category != "surface"):
                    for column_index, plot_type in enumerate(["deltat", "similarity", "cc"]):
                        ax = figs.add_subplot(
                            8, 3, (row_index+1) * 3 + column_index+1)
                        for interation_index, each_iteration in enumerate(iterations_list):
                            sns.distplot(collecction_all[plot_type][interation_index], ax=ax, hist=False, label=f"before iteration {each_iteration}",
                                         kde_kws={"linewidth": 6})
                            if (plot_type == "deltat"):
                                ax.set_xlim((-10, 10))
                            elif(plot_type == "similarity"):
                                ax.set_xlim((0, 1))
                            elif(plot_type == "cc"):
                                ax.set_xlim((0, 1))
                        if (column_index == 0):
                            ax.get_yaxis().set_ticklabels([])
                            ax.set_ylabel(
                                "all phases", fontsize=50, rotation=90)
                        else:
                            ax.get_yaxis().set_ticklabels([])
                        ax.tick_params(axis="x", labelsize=30)
                        if(plot_type != "similarity"):
                            ax.set_xlabel(plot_type, fontsize=30)
                        else:
                            ax.set_xlabel("zero-lag cc", fontsize=30)

                figs.tight_layout()
                pdf.savefig(figs)
                plt.close(fig=figs)


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--misfit_windows_base', required=True, type=str, help="the misfit windows base directory")
    @click.option('--cmts_directory', required=True, type=str, help="the cmt solution directory for all the iterations")
    @click.option('--iterations', required=True, type=str, help="the iterations to plot. eg: 1,2,3")
    @click.option('--pdf_fname', required=True, type=str, help="the pdf path to save")
    @click.option('--snr_threshold', required=True, type=float, help="the snr threshold, such as 3")
    @click.option('--data_info_directory', required=True, type=str, help="the data info directory")
    def main(misfit_windows_base, cmts_directory, iterations, pdf_fname, snr_threshold, data_info_directory):
        iterations_list = iterations.split(",")
        misfit_windows_collection = load_misfit_windows(
            misfit_windows_base, iterations_list)
        _, _, _, event_depth = load_first_arrival_baz_evdp(data_info_directory)
        plot_to_pdf(pdf_fname, cmts_directory, misfit_windows_collection,
                    iterations_list, snr_threshold, event_depth)
    main()
