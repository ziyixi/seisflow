"""
build_structure.py: build the specfem3D structure.
"""
from ...tasks.utils.structure import init_structure


class Build_structure(object):
    def __init__(self, base=None, cmtfiles=None, ref=None, output=None, database=None):
        super().__init__()
        self.base = base
        self.cmtfiles = cmtfiles
        self.ref = ref
        self.output = output
        self.database = database

    def run(self):
        init_structure(self.base, self.cmtfiles,
                       self.ref, self.output, self.database)


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--base', required=True, type=str, help="the base dir to be run.")
    @click.option('--cmtfiles', required=True, type=str, help='cmt files, each named as the id of the event')
    @click.option('--ref', required=True, type=str, help='reference specfem directories')
    @click.option('--output', required=True, type=str, help='directory to place OUTPUT_FILES')
    @click.option('--database', required=True, type=str, help='directory to place DATABASES_MPI')
    def main(base, cmtfiles, ref, output, database):
        run_script = Build_structure(
            base=base, cmtfiles=cmtfiles, ref=ref, output=output, database=database)
        run_script.run()

    main()  # pylint: disable=no-value-for-parameter
