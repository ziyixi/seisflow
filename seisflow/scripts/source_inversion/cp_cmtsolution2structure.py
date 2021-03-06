"""
cp_cmtsolution2structure.py: copy the cmtsolution files in a directory to the simulation directory.
"""
from glob import glob
from os.path import basename, join

import sh

if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--cmtsolution_directory', required=True, type=str, help="the cmtsolution files directory")
    @click.option('--base_directory', required=True, type=str, help="the simulation base directory")
    def main(cmtsolution_directory, base_directory):
        all_files = sorted(glob(join(cmtsolution_directory, "*")))
        for each_file in all_files:
            gcmtid = basename(each_file)
            to_path = join(base_directory, gcmtid, "DATA", "CMTSOLUTION")
            try:
                sh.cp(each_file, to_path)
            except sh.ErrorReturnCode_1:
                pass
    main()  # pylint: disable=no-value-for-parameter
