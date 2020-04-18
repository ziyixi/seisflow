from ...tasks.asdf.seed2asdf import seed2asdf


class Seed2Asdf():
    def __init__(self, seed_directory, cmt_path, output_path):
        super().__init__()
        self.seed_directory = seed_directory
        self.cmt_path = cmt_path
        self.output_path = output_path

    def run(self):
        seed2asdf(self.seed_directory, self.cmt_path, self.output_path)


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--seed_directory', required=True, type=str, help="the directory that contains all the seed files")
    @click.option('--cmt_path', required=True, type=str, help="the cmt solution path")
    @click.option('--output_path', required=True, type=str, help="the output asdf path")
    def main(seed_directory, cmt_path, output_path):
        run_script = Seed2Asdf(seed_directory, cmt_path, output_path)
        run_script.run()

    main()  # pylint: disable=no-value-for-parameter
