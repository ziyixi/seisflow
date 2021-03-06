"""
move_adjoint2sem.py: move the adjoint asdf files back to the SEM folders.
"""
from glob import glob
from os.path import basename, join

import sh


def get_all_adjoint_paths(adjoint_directoy):
    all_files = sorted(glob(join(adjoint_directoy, "*h5")))
    return all_files


def ln_all_files(all_files, base_directory):
    for each_file in all_files:
        gcmtid = basename(each_file).split(".")[0]
        to_path = join(base_directory, gcmtid, "SEM", "adjoint.h5")
        # in the second iteration, the already existing adjoint.h5 will prevent linking
        try:
            sh.unlink(to_path)
        except sh.ErrorReturnCode_1:
            pass
        sh.ln("-s", each_file, to_path)


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--adjoint_directoy', required=True, type=str, help="the directory of adjoint sources")
    @click.option('--base_directory', required=True, type=str, help="the simulation directory")
    def main(adjoint_directoy, base_directory):
        all_files = get_all_adjoint_paths(adjoint_directoy)
        ln_all_files(all_files, base_directory)
    main()  # pylint: disable=no-value-for-parameter
