import argparse
from glob import glob
from os.path import join, basename
import sys

import numpy as np
import sh
import subprocess


def get_args(args=None):
    parser = argparse.ArgumentParser(
        description='A python script to init the structure of the specfem forward simulation',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--base',
                        help='the directory to place all the specefem directories',
                        required=True)

    parser.add_argument('--cmtfiles',
                        help='cmt files, each named as the id of the event',
                        required=True)

    parser.add_argument('--ref',
                        help='reference specfem directories',
                        required=True)

    parser.add_argument('--output',
                        help='directory to place OUTPUT_FILES',
                        required=True)

    parser.add_argument('--database',
                        help='directory to place DATABASES_MPI',
                        required=True)

    results = parser.parse_args(args)
    # return results["base"], results["cmtfiles"], results["ref"], results["output"]
    return results.base, results.cmtfiles, results.ref, results.output, results.database


def init_structure(base, cmtfiles, ref, output, database):
    """
    only copy or ln the necessary files to each simulation directories.
    """
    # base dir
    sh.mkdir("-p", base)
    sh.mkdir("-p", output)
    sh.mkdir("-p", database)

    # some files in the ref
    ref_bin_path = join(ref, "bin")
    # ref_station = join(ref, "DATA", "STATIONS")
    # ref_parfile = join(ref, "DATA", "Par_file")
    # ref_gll = join(ref, "DATA", "GLL")
    ref_values_from_mesher = join(ref, "OUTPUT_FILES", "values_from_mesher.h")

    all_gcmt_ids = glob(join(cmtfiles, "*"))
    for each_gcmtid_path in all_gcmt_ids:
        each_gcmtid = basename(each_gcmtid_path)
        # make running directory
        working_dir = join(base, each_gcmtid)
        sh.mkdir("-p", working_dir)
        # handle DATA
        sh.mkdir("-p", join(working_dir, "DATA"))
        subprocess.call(
            f"ln -s {join(ref, 'DATA', '*')} {join(working_dir, 'DATA/')}", shell=True)
        sh.unlink(join(working_dir, "DATA", "CMTSOLUTION"))
        sh.cp(each_gcmtid_path, join(working_dir, "DATA", "CMTSOLUTION"))
        # handle DATABASES_MPI
        sh.mkdir("-p", join(database, each_gcmtid))
        sh.ln("-s", join(database, each_gcmtid),
              join(working_dir, "DATABASES_MPI"))
        # handle OUTPUT_FILES
        sh.mkdir("-p", join(output, each_gcmtid))
        sh.cp(ref_values_from_mesher, join(
            output, each_gcmtid, "values_from_mesher.h"))
        sh.ln("-s", join(output, each_gcmtid),
              join(working_dir, "OUTPUT_FILES"))
        # handle bin files
        sh.ln("-s", ref_bin_path, join(working_dir, "bin"))


if __name__ == "__main__":
    base, cmtfiles, ref, output, database = get_args(sys.argv[1:])
    init_structure(base, cmtfiles, ref, output, database)
