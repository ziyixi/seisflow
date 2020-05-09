"""
cp_sbatch_to_all.py: cp all the sbatch dirs to the simulation directory.
"""
from glob import glob
from os.path import isdir, join

import click
import sh


@click.command()
@click.option('--source_dir', required=True, type=str, help="the source sbatch dir")
@click.option('--base_dir', required=True, type=str, help="the base simulation dir")
def main(source_dir, base_dir):
    all_paths = sorted(glob(join(base_dir, "*")))
    for each_path in all_paths:
        target = join(each_path, "bin_knl")
        if (not isdir(target)):
            sh.cp("-r", source_dir, target)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
