from ..asdf.sac2asdf import sac2asdf


class Sac2Asdf():
    def __init__(self, sac_directory, response_directory, cmt_path, output_path):
        super().__init__()
        self.sac_directory = sac_directory
        self.response_directory = response_directory
        self.cmt_path = cmt_path
        self.output_path = output_path

    def run(self):
        sac2asdf(self.sac_directory, self.response_directory,
                 self.cmt_path, self.output_path)


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--sac_directory', required=True, type=str, help="the directory that contains all the sac files")
    @click.option('--response_directory', required=True, type=str, help="the directory that contains all the resp files")
    @click.option('--cmt_path', required=True, type=str, help="the cmt solution path")
    @click.option('--output_path', required=True, type=str, help="the output asdf path")
    def main(sac_directory, response_directory, cmt_path, output_path):
        run_script = Sac2Asdf(
            sac_directory, response_directory, cmt_path, output_path)
        run_script.run()

    main()  # pylint: disable=no-value-for-parameter
