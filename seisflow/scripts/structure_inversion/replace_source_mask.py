"""
Replace the source mask.
"""
from glob import glob
from os.path import join, basename
import subprocess
import tqdm
import click


@click.command()
@click.option('--base_directory', required=True, type=str, help="the base simulation directory")
@click.option('--source_mask_directory', required=True, type=str, help="the source mask directory")
def replace_source(base_directory, source_mask_directory):
    all_events_dir = glob(join(base_directory, "*"))
    for each_event_dir in tqdm.tqdm(all_events_dir):
        gcmtid = basename(each_event_dir)
        database_dir = join(each_event_dir, "DATABASES_MPI")
        source_mask_dir = join(source_mask_directory, gcmtid)
        raw_source_masks = join(database_dir, "proc*_reg1_mask_source.bin")
        new_source_masks = join(source_mask_dir, "proc*_reg1_mask_source.bin")
        # rm raw
        subprocess.call(f"rm -rf {raw_source_masks} ", shell=True)
        # ln new
        subprocess.call(f"ln -s {new_source_masks} {database_dir}", shell=True)


if __name__ == "__main__":
    replace_source()  # pylint: disable=no-value-for-parameter
