from .slurmpy import Slurm


def submit_job(jobname, thecommand, N_node, ntasks, partition, time, account, depends_on=None):
    """
    submit_job: a wrapper for slurmpy
        + N_node: total number of nodes used.
        + ntasks: total number of mpi processes.
        + partition: partition used, eg: skx-normal.
        + time: used in slurm format.
        + account: account used in the slurm system.

    """
    if(hpc == "stampede2"):
        s = Slurm(jobname, {"nodes": N_node, "ntasks": ntasks,
                            "partition": partition, "time": time, "account": account})
    elif (hpc == "icer"):
        # virasdf may consume so much memory (estimate as 10*single, around 5G)
        s = Slurm(jobname, {"nodes": N_node, "ntasks": ntasks,
                            "time": time, "cpus-per-task": 1, "mem-per-cpu": "5G"})
    else:
        raise Exception(
            "not supported hpc platform, can be either stampede2 or ICER.")
    job_id = s.run(thecommand, depends_on=depends_on)
    return job_id
