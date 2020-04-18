import subprocess
from glob import glob
from os.path import basename, join

import sh


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

    all_gcmt_ids = sorted(glob(join(cmtfiles, "*")))
    for each_gcmtid_path in all_gcmt_ids:
        each_gcmtid = basename(each_gcmtid_path)
        # make running directory
        working_dir = join(base, each_gcmtid)
        sh.mkdir("-p", working_dir)
        # handle DATA
        sh.mkdir("-p", join(working_dir, "DATA"))
        subprocess.call(
            f"ln -s {join(ref, 'DATA', '*')} {join(working_dir, 'DATA/')}", shell=True)
        sh.unlink(join(working_dir, "DATA", "CMTSOLUTION")
                  )
        sh.unlink(join(working_dir, "DATA", "STATIONS")
                  )
        sh.unlink(join(working_dir, "DATA", "Par_file")
                  )
        try:
            sh.unlink(join(working_dir, "DATA", "STATIONS_ADJOINT"))
        except sh.ErrorReturnCode_1:
            pass
        sh.cp(each_gcmtid_path, join(
            working_dir, "DATA", "CMTSOLUTION"))
        sh.cp(join(ref, 'DATA', 'STATIONS'), join(
            working_dir, "DATA", "STATIONS"))
        sh.cp(join(ref, 'DATA', 'Par_file'), join(
            working_dir, "DATA", "Par_file"))
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
        # make SEM directory
        sh.mkdir("-p", join(working_dir, "SEM"))
        # copy change_simulation_type.pl script for changing the simulation type
        sh.cp(join(ref, "change_simulation_type.pl"), join(
            working_dir, "change_simulation_type.pl"))
        # handle bin files
        sh.ln("-s", ref_bin_path, join(working_dir, "bin"))
