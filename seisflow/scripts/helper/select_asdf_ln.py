"""
select_asdf.py: select asdf files based on the cmt solution files.
"""
from glob import glob
from os.path import join, basename
import subprocess
import click
from os.path import isfile


@click.command()
@click.option('--cmt_directory', required=True, type=str, help="the cmts directory")
@click.option('--asdf_directory', required=True, type=str, help="the source asdf directory")
@click.option('--output_directory', required=True, type=str, help="the destination asdf directory")
def main(cmt_directory, asdf_directory, output_directory):
    all_cmt_files = glob(join(cmt_directory, "*"))
    for each_cmt_file in all_cmt_files:
        gcmtid = basename(each_cmt_file)
        from_path = join(asdf_directory, f"{gcmtid}.h5")
        to_path = join(output_directory, f"{gcmtid}.h5")
        command = f"ln -s {from_path} {to_path} "
        if(isfile(from_path)):
            subprocess.call(command, shell=True)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
