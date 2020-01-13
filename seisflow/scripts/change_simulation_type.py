"""
change_simulation_type.py: change the simulation type in specfem massively.
"""

from glob import glob
from os.path import join, basename
import sh

if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--base_directory', required=True, type=str, help="the simulation directory")
    @click.option('--simulation_type', required=True, type=str, help="can be forward, source, structure and forward_save")
    def main(base_directory, simulation_type):
        flag = None
        if(simulation_type == "forward"):
            flag = "-f"
        elif(simulation_type == "source"):
            flag = "-a"
        elif(simulation_type == "structure"):
            flag = "-b"
        elif(simulation_type == "forward_save"):
            flag = "-F"
        else:
            raise Exception("no such simulation type")
        all_simulation_directories = sorted(glob(join(base_directory, "*")))
        current_path = str(sh.pwd())[:-1]
        for each_simulation_directory in all_simulation_directories:
            sh.cd(each_simulation_directory)
            sh.perl("change_simulation_type.pl", flag)
            sh.cd(current_path)
    main()
