"""
ln_change_tags.py: ln bin files and erase the corresponding tags.
"""
from glob import glob
from os.path import basename, join

import click
import sh


@click.command()
@click.option('--raw_directory', required=True, type=int, help="the directory where bin files are with the appending")
@click.option('--output_directory', required=True, type=int, help="the output directory")
@click.option('--tags', required=True, type=int, help="the bin file tags to select")
def main(raw_directory, output_directory, tags):
    """
    main: ln bin files and erase the corresponding tags.
    """
    tags = tags.split(",")
    for each_tag in tags:
        all_files = sorted(glob(join(raw_directory, f"*_{each_tag}_*")))
        all_basenames = [basename(item) for item in all_files]
        for each_fname in all_basenames:
            each_fname_splitted = each_fname.split("_")
            new_fname = "_".join(each_fname_splitted[:-1]) + ".bin"
            from_path = join(raw_directory, each_fname)
            to_path = join(output_directory, new_fname)
            sh.ln("-s", from_path, to_path)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
