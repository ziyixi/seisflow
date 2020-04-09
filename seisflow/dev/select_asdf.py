"""
select_asdf.py: select asdf files based on the cmt solution files.
"""
from glob import glob
from os.path import join, basename
import subprocess
import click


@click.command()
@click.option('--cmt_directory', required=True, type=str, help="the cmts directory")
@click.option('--asdf_directory', required=True, type=str, help="the source asdf directory")
@click.option('--output_directory', required=True, type=str, help="the destination asdf directory")
def main(cmt_directory, asdf_directory, output_directory):
    all_cmt_files = glob(join(cmt_directory, "*"))
    cp_list = []
    for each_cmt_file in all_cmt_files:
        gcmtid = basename(each_cmt_file)
        if(gcmtid not in cp_list):
            cp_list.append(gcmtid)
            from_path = join(asdf_directory, f"{gcmtid}*")
            to_path = output_directory
            command = f"cp {from_path} {to_path} "
            subprocess.call(command, shell=True)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
