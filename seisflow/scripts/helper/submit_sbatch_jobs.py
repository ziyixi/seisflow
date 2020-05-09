"""
submit_sbatch_jobs.py: submit sbatch jobs in the base simulation directory.
"""
from glob import glob
from os.path import basename, join

import click
import sh


@click.command()
@click.option('--base_dir', required=True, type=str, help="the simulation base dir")
@click.option('--cmd_name', required=True, type=str, help="the cmd name to sbatch")
@click.option('--id_range', required=True, type=str, help="range of events, min,max")
def main(base_dir, cmd_name, id_range):
    all_dirs = sorted(glob(join(base_dir, "*")))
    id_min, id_max = map(int, id_range.split(","))
    used_dirs = all_dirs[id_min:id_max]
    for each_dir in used_dirs:
        sbatch_dir = join(each_dir, "sbatch")
        sh.cd(sbatch_dir)
        sh.sbatch(cmd_name)
        print(f"submitted {cmd_name} in {basename(each_dir)}")


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
