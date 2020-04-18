"""
cp_stations_adjoint2structure.py: copy the STATIONS_ADJOINT files in a directory to the simulation directory.
"""
from glob import glob
from os.path import basename, join

import sh

if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--stations_adjoint_directory', required=True, type=str, help="the stations_adjoint files directory")
    @click.option('--base_directory', required=True, type=str, help="the simulation base directory")
    def main(stations_adjoint_directory, base_directory):
        all_files = sorted(glob(join(stations_adjoint_directory, "*")))
        for each_file in all_files:
            gcmtid = basename(each_file)
            to_path = join(base_directory, gcmtid, "DATA", "STATIONS_ADJOINT")
            sh.cp(each_file, to_path)
    main()  # pylint: disable=no-value-for-parameter
