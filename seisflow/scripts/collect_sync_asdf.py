"""
collect_sync_asdf.py: collect all the synthetics in the output directory.
"""
from os.path import join
from glob import glob
import sh


def get_all_sync_paths(search_directoy):
    all_files = sorted(glob(join(search_directoy, "*", "synthetic.h5")))
    return all_files


def move_all_files(all_files, output_directory):
    sh.mkdir("-p", output_directory)
    for each_file in all_files:
        fname_splitter = each_file.split("/")
        gcmtid = fname_splitter[-2]
        output_path = join(output_directory, f"{gcmtid}.h5")
        sh.mv(each_file, output_path)


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--search_directoy', required=True, type=str, help="the output directory of the specfem result")
    @click.option('--output_directory', required=True, type=str, help="the output directory of the asdf files")
    def main(search_directoy, output_directory):
        all_files = get_all_sync_paths(search_directoy)
        move_all_files(all_files, output_directory)
    main()
