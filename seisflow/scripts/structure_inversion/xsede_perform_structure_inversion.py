"""
xsede_perform_structure_inversion.py: perform the structure inversion for the second iteration and later.
"""
import sys
from os.path import join

import click
import sh

from ..build_structure import Build_structure
from ...tasks import forward_task


def main(base_directory, cmts_directory, ref_directory,
         n_total, n_each, n_iter, nproc):
    """
    main: the main function combining all the parts together.
    """
    # * we have to build the structure to perform the structure inversion.
    build_inversion_structure(base_directory, cmts_directory, ref_directory)
    # * here we have to init the slurm script, no need to load modules here
    result = "date; \n"
    pyexec = sys.executable
    current_path = str(sh.pwd())[:-1]  # pylint: disable=not-callable
    result += f"cd {current_path}; \n"
    # * change the flags to -F
    result += f"{pyexec}  -m seisflow.scripts.change_simulation_type --base_directory {join(base_directory, 'simulation')} --simulation_type forward_save; \n"
    result += f"cd {current_path}; \n"
    # * submit the forward simulation job
    forward_simulation_command = forward_task(base=join(base_directory, "simulation"),
                                              N_total=n_total, N_each=n_each, N_iter=n_iter, nproc=nproc, run_mesh=True)
    result += forward_simulation_command
    result += f"cd {current_path}; \n"
    # * collect the sync from the forward simulation
    result += f"{pyexec} -m seisflow.scripts.collect_sync_asdf --search_directoy {join(base_directory, 'output')} --output_directory {join(base_directory, 'raw_sync')}; \n"
    result += f"cd {current_path}; \n"


def build_inversion_structure(base_directory, cmts_directory, ref_directory):
    """
    build_inversion_structure: build the structure to contain all the essencial directories used in the inversion and the simulation directory.
    """
    sh.mkdir("-p", base_directory)
    # * copy cmts_directory
    sh.cp("-r", cmts_directory, join(base_directory, "cmts"))
    # * init the simulation directory
    output_path = join(base_directory, "output")
    sh.mkdir("-p", output_path)
    database_path = join(base_directory, "database")
    sh.mkdir("-p", database_path)
    simulation_path = join(base_directory, "simulation")
    sh.mkdir("-p", simulation_path)
    run_script = Build_structure(
        base=simulation_path, cmtfiles=join(base_directory, "cmts"), ref=ref_directory,
        output=output_path, database=database_path)
    run_script.run()
    # * make the directory for the sync of the forward simulation
    sh.mkdir("-p", join(base_directory, "raw_sync"))
