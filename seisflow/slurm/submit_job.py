from .slurmpy import Slurm


def submit_job(thecommand, N_node, ntasks, partition, time, account):
    """
    submit_job: a wrapper for slurmpy
        + N_node: total number of nodes used.
        + ntasks: total number of mpi processes.
        + partition: partition used, eg: skx-normal.
        + time: used in slurm format.
        + account: account used in the slurm system.

    """
    s = Slurm("sync", {"nodes": N_node, "ntasks": ntasks,
                       "partition": partition, "time": time, "account": account})
    job_id = s.run(thecommand)
    return job_id
