import argparse
import sys
from glob import glob
from os.path import join, dirname, abspath

import configparser


def get_args(args=None):
    parser = argparse.ArgumentParser(
        description='A python script to submit jobs in one sbatch job',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--base',
                        help='the directory to place all the specefem directories',
                        required=True)
    parser.add_argument('--configure',
                        help='the configuration file name in the configuration directory',
                        required=True)
    results = parser.parse_args(args)
    current_file_dir = dirname(abspath(__file__))
    return results.base, join(current_file_dir, "..", "configuration", results.configure)


def get_dirs(base):
    return glob(join(base, "*"))


def get_scripts(thedirs, N_total, N_each, N_iter, nproc):
    result = "date; "
    result += "module load boost/1.68; "
    result += "module load phdf5/1.8.16; "
    # for xmeshfem3D
    result += f"echo 'start xmeshfem3D'; "
    for iiter in range(N_iter):
        result += f"echo 'start iteration {iiter}'; "
        for ieach in range(N_each):
            # ievent
            ievent = iiter * N_each + ieach
            if (ievent >= N_total):
                continue
            ievent_dir = thedirs[ievent]
            # cd
            result += f"cd {ievent_dir}; "
            # if N_each-1
            if(ieach == N_each-1):
                inc = ieach*nproc
                result += f"ibrun -n {nproc} -o {inc} ./bin/xmeshfem3D; "
            else:
                inc = ieach*nproc
                result += f"ibrun -n {nproc} -o {inc} ./bin/xmeshfem3D & "
        result += f"wait; "
        result += f"echo 'end iteration {iiter}'; "
        result += f"date; "

    # for xspecfem3D
    result += f"echo 'start xspecfem3D'; "
    for iiter in range(N_iter):
        result += f"echo 'start iteration {iiter}'; "
        for ieach in range(N_each):
            # ievent
            ievent = iiter * N_each + ieach
            if (ievent >= N_total):
                continue
            ievent_dir = thedirs[ievent]
            # cd
            result += f"cd {ievent_dir}; "
            # if N_each-1
            if(ieach == N_each-1):
                inc = ieach*nproc
                result += f"ibrun -n {nproc} -o {inc} ./bin/xspecfem3D; "
            else:
                inc = ieach*nproc
                result += f"ibrun -n {nproc} -o {inc} ./bin/xspecfem3D & "
        result += f"wait; "
        result += f"echo 'end iteration {iiter}'; "
        result += f"date; "

    return result


def forward_task(base=None, N_total=None, N_each=None, N_iter=None, nproc=None):
    """
    submit a forward job of specfem3d-globe.
        + base: the base dir to be run.
        + N_total: total number of events.
        + N_each: number of running jobs at each iterations.
        + N_iter: number of iterations to run.
        + nproc: number of mpi processes for the each event.
    """
    thedirs = get_dirs(base)
    thecommand = get_scripts(thedirs, N_total, N_each, N_iter,
                             nproc)
    return thecommand
