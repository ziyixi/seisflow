"""
sac2asdf_hinet.py: convert the sac files downloaded from hinetpy to the asdf format.
"""
from ...tasks.asdf.sac2asdf_hinet import sac2asdf_hinet


class Sac2Asdf():
    def __init__(self, sac_directory, cmt_path, output_path):
        super().__init__()
        self.sac_directory = sac_directory
        self.cmt_path = cmt_path
        self.output_path = output_path

    def run(self):
        sac2asdf_hinet(self.sac_directory, self.cmt_path, self.output_path)


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--sac_directory', required=True, type=str, help="the directory that contains all the sac files")
    @click.option('--cmt_path', required=True, type=str, help="the cmt solution path")
    @click.option('--output_path', required=True, type=str, help="the output asdf path")
    def main(sac_directory,  cmt_path, output_path):
        run_script = Sac2Asdf(
            sac_directory, cmt_path, output_path)
        run_script.run()

    main()  # pylint: disable=no-value-for-parameter
