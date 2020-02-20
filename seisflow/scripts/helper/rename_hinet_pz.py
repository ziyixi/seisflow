"""
rename_hinet_pz.py: rename hinet pz files in order to process.
"""
from os.path import join, basename, dirname
from glob import glob
import sh


def rename_pz_hinet(pz_directory):
    """
    rename pz files in a directory. (to net.sta.cha)
    """
    all_files = sorted(glob(join(pz_directory, "*")))
    for each_file in all_files:
        fname = basename(each_file)
        net, sta, cha, _ = fname.split(".")
        cha = "HH" + cha
        if (cha == "HHU"):
            cha = "HHZ"
        newfname = f"{net}.{sta}.{cha}"
        new_file = join(dirname(each_file), newfname)
        sh.mv(each_file, new_file)


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--pz_directory', required=True, type=str, help="the pz directory.")
    def main(pz_directory):
        """
        rename hinet pz files in order to process.
        """
        rename_pz_hinet(pz_directory)

    main()  # pylint: disable=no-value-for-parameter
