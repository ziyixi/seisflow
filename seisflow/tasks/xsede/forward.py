from glob import glob
from os.path import join


def get_dirs(base):
    return sorted(glob(join(base, "*")))


def get_scripts(thedirs, N_total, N_each, N_iter, nproc, run_mesh=True):
    result = "date; "
    result += "module load boost/1.68; "
    result += "module load phdf5/1.8.16; "
    # for xmeshfem3D
    if(run_mesh):
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
            result += f"date; \n"

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
        result += f"date; \n"

    return result


def forward_task(base=None, N_total=None, N_each=None, N_iter=None, nproc=None, run_mesh=True):
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
                             nproc, run_mesh=run_mesh)
    return thecommand
