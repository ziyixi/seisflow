"""
Calculate the overall misfit from the misfits directory and the weight file. 
"""
from glob import glob
from os.path import basename, join

import click
from seisflow.utils.load_files import load_pickle


def read_files(misfit_directory,weight_file):
    weights=load_pickle(weight_file)
    misfits={}
    misfit_files=sorted(glob(join(misfit_directory,"*pkl")))
    for each_file in misfit_files:
        thekey=basename(each_file).split(".")[0]
        misfits[thekey]=load_pickle(each_file)
    return weights,misfits

def cal_overall_misfit(weights,misfits):
    overall_misfit=0
    overall_weight=0
    # 6 categories
    misfit_categories={
        "z":0,
        "r":0,
        "t":0,
        "surface_z":0,
        "surface_r":0,
        "surface_t":0
    }
    weight_categories={
        "z":0,
        "r":0,
        "t":0,
        "surface_z":0,
        "surface_r":0,
        "surface_t":0
    }

    for each_event in weights:
        for each_station in weights[each_event]:
            for each_category in weights[each_event][each_station]:
                for iwin,each_window in enumerate(weights[each_event][each_station][each_category]):
                    theweight=each_window.snr*each_window.cc*each_window.deltat*each_window.geographical*each_window.category
                    misfit_win=misfits[each_event][each_station][each_category].windows[iwin]
                    themisfit=1-misfit_win.similarity
                    overall_misfit+=theweight*themisfit
                    overall_weight+=theweight
                    misfit_categories[each_category]+=theweight*themisfit
                    weight_categories[each_category]+=theweight
    overall_misfit=overall_misfit/overall_weight
    for each_category in ["z","r","t","surface_z","surface_r","surface_t"]:
        misfit_categories[each_category]=misfit_categories[each_category]/weight_categories[each_category]
    return overall_misfit,misfit_categories["z"],misfit_categories["r"],misfit_categories["t"],misfit_categories["surface_z"],misfit_categories["surface_r"],misfit_categories["surface_t"]

@click.command()
@click.option('--misfit_directory', required=True, type=str, help="the misfit directory")
@click.option('--weight_file', required=True, type=str, help="the weight file")
def main(misfit_directory,weight_file):
    weights,misfits=read_files(misfit_directory,weight_file)
    result=cal_overall_misfit(weights,misfits)
    print(result)

if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
    