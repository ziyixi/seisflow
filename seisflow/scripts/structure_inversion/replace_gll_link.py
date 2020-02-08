"""
replace_gll_link.py: replace all the links in the forward simulation directory.
"""
from os.path import join
from glob import glob
import sh


def relink_single(single_event_directory, new_gll_directory):
    """
    relink GLL for each simulation directory
    """
    gll_path = join(single_event_directory, "DATA", "GLL")
    # * definitely there will be a GLL
    sh.unlink(gll_path)
    sh.ln("-s", new_gll_directory, gll_path)


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--simulation_directory', required=True, type=str, help="the simulation directory")
    @click.option('--new_gll_directory', required=True, type=str, help="the new gll directory")
    def main(simulation_directory, new_gll_directory):
        """
        relink GLL for all events simulation directory
        """
        all_events_directories = sorted(glob(join(simulation_directory, "*")))
        for each_event_directory in all_events_directories:
            relink_single(each_event_directory, new_gll_directory)
    main()  # pylint: disable=no-value-for-parameter
