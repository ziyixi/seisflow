"""
A scripts to build the forward simulation structure and submit the job.
"""
from ...slurm.submit_job import submit_job
from ...tasks.xsede.forward import forward_task

# * test is passed for this script on 01/07/2020


class Run_multiple_forward_jobs(object):
    def __init__(self, base=None, N_total=None, N_each=None, N_iter=None,
                 nproc=None, N_node=None, ntasks=None, partition=None, time=None, account=None, run_mesh=True):
        super().__init__()
        self.base = base
        self.N_total = N_total
        self.N_each = N_each
        self.N_iter = N_iter
        self.nproc = nproc
        self.N_node = N_node
        self.ntasks = ntasks
        self.partition = partition
        self.time = time
        self.account = account
        self.run_mesh = run_mesh

    def run(self):
        thecommand = forward_task(
            base=self.base, N_total=self.N_total, N_each=self.N_each, N_iter=self.N_iter, nproc=self.nproc, run_mesh=self.run_mesh)
        job_id = submit_job("forward", thecommand, self.N_node, self.ntasks,
                            self.partition, self.time, self.account, "stampede2")
        return job_id


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--base', required=True, type=str, help="the base dir to be run.")
    @click.option('--ntotal', required=True, type=int, help="total number of events.")
    @click.option('--neach', required=True, type=int, help="number of running jobs at each iterations.")
    @click.option('--niter', required=True, type=int, help="number of iterations to run.")
    @click.option('--nproc', required=True, type=int, help="number of mpi processes for each event.")
    @click.option('--nnode', required=True, type=int, help="total number of nodes used.")
    @click.option('--ntasks', required=True, type=int, help="total number of mpi processes.")
    @click.option('--partition', required=True, type=str, help="partition used, eg: skx-normal.")
    @click.option('--time', required=True, type=str, help="used in slurm format.")
    @click.option('--account', required=True, type=str, help="account used in the slurm system.")
    @click.option('--run_mesh/--no-run_mesh', required=True, default=True)
    def main(base, ntotal, neach, niter, nproc, nnode, ntasks, partition, time, account, run_mesh):
        run_script = Run_multiple_forward_jobs(base=base,  N_total=ntotal,
                                               N_each=neach, N_iter=niter, nproc=nproc, N_node=nnode, ntasks=ntasks, partition=partition, time=time, account=account, run_mesh=run_mesh)
        run_script.run()

    main()  # pylint: disable=no-value-for-parameter
