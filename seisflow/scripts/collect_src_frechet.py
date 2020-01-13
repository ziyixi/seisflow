"""
collect_src_frechet.py: collect all the src_frechet files in the output directory.
"""
from glob import glob
from os.path import join

import sh


def get_all_src_frechet_paths(search_directoy):
    all_files = sorted(glob(join(search_directoy, "*", "src_frechet.000001")))
    return all_files


def move_all_files(all_files, output_directory):
    sh.mkdir("-p", output_directory)
    for each_file in all_files:
        fname_splitter = each_file.split("/")
        gcmtid = fname_splitter[-2]
        output_path = join(output_directory, f"{gcmtid}")
        sh.mv(each_file, output_path)


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--search_directoy', required=True, type=str, help="the output directory of the specfem result")
    @click.option('--output_directory', required=True, type=str, help="the output directory of the src_frechet files")
    def main(search_directoy, output_directory):
        all_files = get_all_src_frechet_paths(search_directoy)
        move_all_files(all_files, output_directory)
    main()
