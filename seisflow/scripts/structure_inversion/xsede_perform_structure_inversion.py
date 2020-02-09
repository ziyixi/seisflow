"""
xsede_perform_structure_inversion.py: perform the structure inversion for the second iteration and later.
"""
import sys
from os.path import join

import click
import sh

from ...slurm.submit_job import submit_job
from ...tasks import forward_task
from ...utils.setting import LINE_SEARCH_PERTURBATION
from ..build_structure import Build_structure
from ..xsede_perform_source_inversion import (calculate_misfit_windows,
                                              calculate_stations_adjoint,
                                              change_simulation_type,
                                              collect_sync_files,
                                              cp_stations_adjoint2structure,
                                              ln_adjoint_source_to_structure)
from ..xsede_process_sync import kernel as get_process_sync_scripts
from .process_kernel import \
    construct_structure as construct_process_kernel_structure
from .process_kernel import kernel as process_kernel


@click.command()
@click.option('--base_directory', required=True, type=str, help="the base inversion directory")
@click.option('--cmts_directory', required=True, type=str, help="the cmts directory")
@click.option('--ref_directory', required=True, type=str, help="the reference specfem directory")
@click.option('--windows_directory', required=True, type=str, help="the windows directory")
@click.option('--data_asdf_directory', required=True, type=str, help="the processed data directory")
@click.option('--data_info_directory', required=True, type=str, help="the data info directory")
@click.option('--last_step_kernel_directory', required=True, type=str, help="the last step smoothed kernel directory")
@click.option('--stations_path', required=True, type=str, help="the stations path")
@click.option('--sem_utils_directory', required=True, type=str, help="the sem_utils directory")
@click.option('--past_raw_directory', required=True, type=str, help="raw directory in the last step")
@click.option('--n_total', required=True, type=int, help="the total number of events")
@click.option('--n_each', required=True, type=int, help="number of events to run in each iteration")
@click.option('--n_iter', required=True, type=int, help="the number of iterations to run")
@click.option('--nproc', required=True, type=int, help="the number of processes used for each event")
@click.option('--n_node', required=True, type=int, help="the number of nodes used in simulation")
@click.option('--partition', required=True, type=str, help="the partion name, eg: skx-normal")
@click.option('--time_forward', required=True, type=str, help="the time used in step 1")
@click.option('--account', required=True, type=str, help="the stampede2 account")
@click.option('--n_node_process_kernel', required=True, type=int, help="the node number in step 2")
@click.option('--time_process_kernel', required=True, type=str, help="the time used in step 2")
@click.option('--time_run_perturbation', required=True, type=str, help="the time used in step 3")
@click.option('--periods', required=True, type=str, help="periods in filtering: minp1,maxp1/minp2,maxp2/...")
@click.option('--waveform_length', required=True, type=int, help="the length of the waveform to cut")
@click.option('--sampling_rate', required=True, type=int, help="the sampling rate to use")
@click.option('--taper_tmin_tmaxs', required=True, type=str, help="the taper time bands: minp1,maxp1/minp2,maxp2/...")
@click.option('--sigma_h', required=True, type=float, help="the value of sigma_h (km)")
@click.option('--sigma_v', required=True, type=float, help="the value of sigma_v (km)")
def main(base_directory, cmts_directory, ref_directory, windows_directory, data_asdf_directory, data_info_directory, last_step_kernel_directory,
         stations_path, sem_utils_directory, past_raw_directory,
         n_total, n_each, n_iter, nproc, n_node, partition, time_forward, account, n_node_process_kernel, time_process_kernel, time_run_perturbation,
         periods, waveform_length, sampling_rate, taper_tmin_tmaxs,
         sigma_h, sigma_v):
    """
    perform the structure inversion for the second iteration and later.
    """
    time = time_forward
    # * we have to build the structure to perform the structure inversion.
    build_inversion_structure(base_directory, cmts_directory, ref_directory)
    # * ======================================================================================================================
    # * here we have to init the slurm script, no need to load modules here
    result = "date; \n"
    pyexec = sys.executable
    current_path = str(sh.pwd())[:-1]  # pylint: disable=not-callable
    # * change the flags to -F
    result += change_simulation_type(pyexec,
                                     join(base_directory, 'simulation'), "forward_save")
    # * submit the forward simulation job
    forward_simulation_command = forward_task(base=join(base_directory, "simulation"),
                                              N_total=n_total, N_each=n_each, N_iter=n_iter, nproc=nproc, run_mesh=True)
    result += forward_simulation_command
    result += f"cd {current_path}; \n"
    # * collect the sync from the forward simulation
    result += collect_sync_files(
        pyexec, {join(base_directory, 'output')}, join(base_directory, 'raw_sync'))
    # * process the sync
    result += get_process_sync_scripts(join(base_directory, "raw_sync"), join(base_directory, "processed_sync"), 1, n_total, nproc,
                                       periods, waveform_length, sampling_rate, taper_tmin_tmaxs, reference_directory=past_raw_directory)
    result += f"cd {current_path}; \n"
    # * calculate the misfit windows
    body_periods, surface_periods = periods.split("/")
    body_periods_splitter = body_periods.split(",")
    surface_periods_splitter = surface_periods.split(",")
    min_periods = f"{body_periods_splitter[0]},{surface_periods_splitter[0]}"
    max_periods = f"{body_periods_splitter[1]},{surface_periods_splitter[1]}"
    result += calculate_misfit_windows(pyexec, n_total,
                                       windows_directory, join(
                                           base_directory, "misfit_windows"), min_periods, max_periods,
                                       data_asdf_directory, join(base_directory, "processed_sync"), data_info_directory)
    # * calculate the adjoint source, and ln it to the sem directory
    result += calculate_adjoint_source(pyexec, n_total,
                                       join(base_directory, "misfit_windows"), stations_path, join(
                                           base_directory, "raw_sync"),
                                       join(
                                           base_directory, "processed_sync"), data_asdf_directory, data_info_directory,
                                       join(base_directory, "adjoint_source"), body_periods, surface_periods)
    result += ln_adjoint_source_to_structure(pyexec,
                                             join(base_directory, "adjoint_source"), join(base_directory, "simulation"))
    # * generate STATIONS_ADJOINT and cp it to the simulation directory
    result += calculate_stations_adjoint(pyexec, stations_path,
                                         join(base_directory, "misfit_windows"), join(base_directory, "stations_adjoint"))
    result += cp_stations_adjoint2structure(pyexec,
                                            join(base_directory, "stations_adjoint"), join(base_directory, "simulation"))
    # * change the simulation type to the type 3
    result += change_simulation_type(pyexec,
                                     join(base_directory, 'simulation'), "structure")
    # * do the adjoint simulation
    adjoint_simulation_command = forward_task(base=join(base_directory, "simulation"),
                                              N_total=n_total, N_each=n_each, N_iter=n_iter, nproc=nproc, run_mesh=False)
    result += adjoint_simulation_command
    result += f"cd {current_path}; \n"
    # * here we submit the first job
    step1_jobid = submit_job("step1_structure", result, n_node, n_each *
                             nproc, partition, time, account, "stampede2")
    # * ======================================================================================================================
    # * now we do the kernel processing part
    # * firstly we have to make the structure (note the async here)
    result = ""
    kernel_process_directory = join(base_directory, "process_kernel")
    input_model_directory = join(ref_directory, "DATA", "GLL")
    construct_process_kernel_structure(
        join(base_directory,
             "database"), ref_directory, sem_utils_directory, kernel_process_directory,
        input_model_directory, last_step_kernel=last_step_kernel_directory)
    # * now perform the structure inversion
    result += process_kernel(join(base_directory,
                                  "process_kernel"), sigma_h, sigma_v, n_total, itern=True)
    # * now submit the job to process the kernel
    step2_jobid = submit_job("step2_structure", result, n_node_process_kernel, n_total, partition, time_process_kernel,
                             account, "stampede2", depends_on=[step1_jobid])
    # * ======================================================================================================================
    # * now we use the generated perturbed model to do a forward simulation
    result = ""
    # * change the flags to -F
    result += change_simulation_type(pyexec,
                                     join(base_directory, 'simulation'), "forward_save")
    result += f"cd {current_path}; \n"
    # * for all events in the simulation dir, change the GLL soft link
    result += replace_gll_link(pyexec, join(base_directory,
                                            "simulation"), join(base_directory, "process_kernel", f"perturbed_{LINE_SEARCH_PERTURBATION}_for_line_search"))
    # * run the forward simulation job
    perturbation_simulation_command = forward_task(base=join(base_directory, "simulation"),
                                                   N_total=n_total, N_each=n_each, N_iter=n_iter, nproc=nproc, run_mesh=True)
    result += perturbation_simulation_command
    result += f"cd {current_path}; \n"
    # * collect the sync in the perturbation directory
    result += collect_sync_files(
        pyexec, {join(base_directory, 'output')}, join(base_directory, 'perturbed_sync'))
    # * process the sync
    result += get_process_sync_scripts(join(base_directory, "perturbed_sync"), join(base_directory, "processed_perturbed_sync"), 1, n_total, nproc,
                                       periods, waveform_length, sampling_rate, taper_tmin_tmaxs)
    result += f"cd {current_path}; \n"
    # * submit the job
    submit_job("step3_structure", result, n_node, n_each *
               nproc, partition, time_run_perturbation, account, "stampede2", depends_on=[step2_jobid])


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
    sh.mkdir("-p", join(base_directory, "processed_sync"))
    # * mkdir for misfit windows
    sh.mkdir("-p", join(base_directory, "misfit_windows"))
    # * mkdir for adjoint source
    sh.mkdir("-p", join(base_directory, "adjoint_source"))
    sh.mkdir("-p", join(base_directory, "stations_adjoint"))
    # * mkdir for kernel processing
    sh.mkdir("-p", join(base_directory, "process_kernel"))
    # * mkdir to collect the perturbed sync
    sh.mkdir("-p", join(base_directory, "perturbed_sync"))
    sh.mkdir("-p", join(base_directory, "processed_perturbed_sync"))


def calculate_adjoint_source(py, nproc, misfit_windows_directory, stations_path, raw_sync_directory, sync_directory,
                             data_directory, data_info_directory, output_directory, body_band, surface_band):
    """
    Calculatet the adjoint source for the structure inversion.
    """
    script = f"ibrun -n {nproc} {py} -m seisflow.scripts.mpi_calculate_adjoint_source_zerolagcc_multiple_events --misfit_windows_directory {misfit_windows_directory} --stations_path {stations_path} --raw_sync_directory {raw_sync_directory} --sync_directory {sync_directory} --data_directory {data_directory} --data_info_directory {data_info_directory} --output_directory {output_directory} --body_band {body_band} --surface_band {surface_band}; \n"
    return script


def replace_gll_link(py, simulation_directory, new_gll_directory):
    """
    replace all gll links.
    """
    script = f"{py} -m seisflow.scripts.structure_inversion.replace_gll_link --simulation_directory {simulation_directory} --new_gll_directory {new_gll_directory}; \n"
    return script


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
