"""
make_green_cmtsolution.py: make the cmtsolution file with half duration as 0.
"""
from glob import glob
from os.path import basename, join

import obspy

from ..utils.save_files import save_cmtsolution

if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--cmtsolution_directory', required=True, type=str, help="the raw cmtsolution directory")
    @click.option('--output_directory', required=True, type=str, help="the green function cmtsolution directory")
    def main(cmtsolution_directory, output_directory):
        """
        make all cmtsolution files in cmtsolution_directory's half duration as 0.
        """
        all_files = glob(join(cmtsolution_directory, "*"))
        for each_file in all_files:
            gcmtid = basename(each_file)
            cmtsolution = obspy.read_events(each_file)[0]
            # use 0 will cause an error in obspy, specfem will change all cmtsolution less that 5*dt as 5*dt.
            cmtsolution.focal_mechanisms[0].moment_tensor.source_time_function.duration = 0.0002
            output_path = join(output_directory, gcmtid)
            save_cmtsolution(output_path, cmtsolution)
    main()
