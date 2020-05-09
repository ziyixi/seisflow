"""
cp_sbatch_to_all.py: cp all the sbatch dirs to the simulation directory.
"""
from os.path import join, isdir
from glob import glob
import sh
import click


@click.command()
@click.option('--source_dir', required=True, type=str, help="the source sbatch dir")
@click.option('--base_dir', required=True, type=str, help="the base simulation dir")
def main(source_dir, base_dir):
    all_paths = sorted(glob(join(base_dir, "*")))
    for each_path in all_paths:
        target = join(each_path, "sbatch")
        if (not isdir(target)):
            sh.cp(source_dir, target)


if __name__ == "__main__":
    main()
