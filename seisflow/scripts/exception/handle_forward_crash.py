"""
handle_forward_crash.py: handle the case when the forward simulation job has crashed due to system's error.
"""
from glob import glob
from os.path import basename, join

import click
import sh


@click.command()
@click.option('--base_directory', required=True, type=str, help="the raw simulation directory")
@click.option('--output_directory', required=True, type=str, help="the new simulation directory")
def main(base_directory, output_directory):
    """
    ln subdirectories that have no synthetic.h5 to a new base directory.
    """
    all_dirs = sorted(glob(join(base_directory, "*")))
    all_syncs = sorted(
        glob(join(base_directory, "*", "OUTPUT_FILES", "synthetic.h5")))
    all_dirs_with_sync = ["/".join(item.split("/")[:-2]) for item in all_syncs]
    all_dirs_no_syncs = sorted(set(all_dirs)-set(all_dirs_with_sync))
    for each_dir in all_dirs_no_syncs:
        thebasename = basename(each_dir)
        from_path = each_dir
        to_path = join(output_directory, thebasename)
        sh.ln("-s", from_path, to_path)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
