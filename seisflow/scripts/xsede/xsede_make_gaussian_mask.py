import click

from ...slurm.submit_job import submit_job
from ..structure_inversion.make_gaussian_mask import make_gaussian_mask


@click.command()
@click.option('--base_dir', required=True, type=str, help="the simulation directory")
@click.option('--sem_bin_dir', required=True, type=str, help="the sem program directory")
@click.option('--output_dir', required=True, type=str, help="the output directory")
@click.option('--nproc', required=True, type=int, help="nproc in the database directory")
@click.option('--n_node', required=True, type=int, help="number of nodes to use")
@click.option('--n_cores', required=True, type=int, help="number of cores to use")
@click.option('--partition', required=True, type=str, help="the partition to use, eg: skx-normal")
@click.option('--used_time', required=True, type=str, help="the running time")
@click.option('--account', required=True, type=str, help="the account to use")
@click.option('--source_mask_km', required=True, type=float, help="the source mask radius in km")
@click.option('--receiver_mask_km', required=True, type=float, help="the receiver mask radius in km")
def main(base_dir, sem_bin_dir, output_dir, nproc, n_node, n_cores, partition, used_time, account, source_mask_km, receiver_mask_km):
    # * firstly we generate mask_xyz files
    thecommand = ""
    thecommand += f"python -m seisflow.scripts.structure_inversion.generate_source_mask_list --base_directory {base_dir} --source_mask_km {source_mask_km} --receiver_mask_km {receiver_mask_km} "
    thecommand += make_gaussian_mask(base_dir, sem_bin_dir, output_dir, nproc)
    submit_job("generate mask", thecommand,
               n_node, n_cores, partition, used_time, account, "stampede2")


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
