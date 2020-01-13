"""
xsede_perform_source_inversion.py: perform source inversion for several events.
"""
from glob import glob
from os.path import join, isdir
import sh
from ..tasks.forward import forward_task
from ..slurm.submit_job import submit_job
import click
from ..tasks.structure import init_structure


def calculate_adjoint_source_raw(py, nproc, misfit_windows_directory, stations_path, raw_sync_directory, sync_directory,
                                 data_directory, data_info_directory, output_directory, body_band, surface_band):
    """
    At the first step, we should calculate the adjoint source for all the events.
    """
    script = f"ibrun -n {nproc} {py} -m seisflow.scripts.mpi_calculate_adjoint_source_zerolagcc --misfit_windows_directory {misfit_windows_directory} --stations_path {stations_path} --raw_sync_directory {raw_sync_directory} --sync_directory {sync_directory} --data_directory {data_directory} --data_info_directory {data_info_directory} --output_directory {output_directory} --body_band {body_band} --surface_band {surface_band}; \n"
    return script


def calculate_stations_adjoint(py, stations_path, misfit_windows_directory, output_directory):
    """
    get the files STATIONS_ADJOINT.
    """
    script = f"ibrun -n 1 {py} -m seisflow.scripts.get_stations_adjoint --stations_path {stations_path} --misfit_windows_directory {misfit_windows_directory} --output_directory {output_directory}; \n"
    return script


def cp_stations_adjoint2structure(py, stations_adjoint_directory, base_directory):
    """
    copy the stations adjoint file to the structure.
    """
    script = f"ibrun -n 1 {py} -m seisflow.scripts.cp_stations_adjoint2structure --stations_adjoint_directory {stations_adjoint_directory} --base_directory {base_directory}; \n"
    return script


def calculate_misfit_windows(py, nproc, windows_directory, output_directory, min_periods, max_periods, data_asdf_directory, sync_asdf_directory, data_info_directory):
    """
    calculate misfit windows
    """
    script = f"ibrun -n {nproc} {py} -m seisflow.scripts.mpi_calculate_misfit_windows --windows_directory {windows_directory} --output_directory {output_directory} --min_periods {min_periods} --max_periods {max_periods} --data_asdf_directory {data_asdf_directory} --sync_asdf_directory {sync_asdf_directory} --data_info_directory {data_info_directory}; \n"
    return script


def ln_adjoint_source_to_structure(py, adjoint_directoy, base_directory):
    """
    cp the adjoint sources to the SEM folders.
    """
    script = f"ibrun -n 1 {py} -m seisflow.scripts.ln_adjoint2sem --adjoint_directoy {adjoint_directoy} --base_directory {base_directory}; \n"
    return script


def generate_green_cmtsolutions(py, cmtsolution_directory, output_directory):
    """
    convert the raw cmtsolution files to tau=0.
    """
    script = f"ibrun -n 1 {py} -m seisflow.scripts.make_green_cmtsolution --cmtsolution_directory {cmtsolution_directory} --output_directory {output_directory}; \n"
    return script


def collect_sync_files(py, search_directoy, output_directory):
    """
    collect sync files from the output directory to the output directory.
    """
    script = f"ibrun -n 1 {py} -m seisflow.scripts.collect_sync_asdf --search_directoy {search_directoy} --output_directory {output_directory}; \n"
    return script


def collect_src_frechet_files(py, search_directoy, output_directory):
    """
    collect src_frechet files
    """
    script = f"ibrun -n 1 {py} -m seisflow.scripts.collect_src_frechet --search_directoy {search_directoy} --output_directory {output_directory}; \n"
    return script


def convert_green2conv_processed(py, nproc, cmts_directory, green_directory, conv_directory, waveform_length, taper_tmin_tmaxs, periods, sampling_rate):
    """
    conv and process the green function.
    """
    script = f"ibrun -n {nproc} {py} -m seisflow.scripts.mpi_convert_green2sync --cmts_directory {cmts_directory} --green_directory {green_directory} --conv_directory {conv_directory} --waveform_length {waveform_length} --taper_tmin_tmaxs {taper_tmin_tmaxs} --periods {periods} --sampling_rate {sampling_rate}; \n"
    return script


def change_simulation_type(py, base_directory, simulation_type):
    """
    change simulation type
    """
    script = f"ibrun -n 1 {py} -m seisflow.scripts.change_simulation_type --base_directory {base_directory} --simulation_type {simulation_type}; \n"
    return script


def make_perturbed_cmtsolution(py, src_frechet_directory, cmtsolution_directory, output_directory):
    """
    make the pertured cmtsolution based on src_frechet.
    """
    script = f"ibrun -n 1 {py} -m seisflow.scripts.make_perturbed_cmtsolution --src_frechet_directory {src_frechet_directory} --cmtsolution_directory {cmtsolution_directory} --output_directory {output_directory}; \n"
    return script


def cp_cmtsolution2structure(py, cmtsolution_directory, base_directory):
    """
    cp cmtsolution files in cmtsolution_directory to the simulation structure.
    """
    script = f"ibrun -n 1 {py} -m seisflow.scripts.cp_cmtsolution2structure --cmtsolution_directory {cmtsolution_directory} --base_directory {base_directory}; \n"
    return script


def source_line_search(py, nproc, green_raw_asdf_directory, green_perturbed_asdf_directory, data_asdf_directory, windows_directory, data_info_directory, stations_path,
                       src_frechet_directory, cmtsolution_directory, output_newcmtsolution_directory,
                       body_band, surface_band, taper_tmin_tmax,
                       alpha_range, t0_range, tau_range):
    script = f"ibrun -n {nproc} {py} -m seisflow.scripts.mpi_source_line_search --green_raw_asdf_directory {green_raw_asdf_directory} --green_perturbed_asdf_directory {green_perturbed_asdf_directory} --data_asdf_directory {data_asdf_directory} --windows_directory {windows_directory} --data_info_directory {data_info_directory} --stations_path {stations_path} --src_frechet_directory {src_frechet_directory} --cmtsolution_directory {cmtsolution_directory} --output_newcmtsolution_directory {output_newcmtsolution_directory} --body_band {body_band} --surface_band {surface_band} --taper_tmin_tmax {taper_tmin_tmax} --alpha_range {alpha_range} --t0_range {t0_range} --tau_range {tau_range}; \n"
    return script


def make_source_inversion_directory(iter_number, inversion_directory, cmtfiles_directory, ref_directory, data_directory, windows_directory,
                                    data_info_directory, stations_path, raw_sync_directory):
    """
    make the running directory for the source inversion.
    """
    # make the structure directory for specfem
    specfem_base = join(inversion_directory, "simulation")
    specfem_cmtfiles = join(inversion_directory, "cmts",
                            "cmtsolutions_current_step")
    specfem_ref = join(inversion_directory, "ref")
    specfem_output = join(inversion_directory, "output")
    specfem_database = join(inversion_directory, "database")
    specfem_stations_adjoint = join(inversion_directory, "stations_adjoint")
    data_processed = join(inversion_directory, "data_processed")
    windows = join(inversion_directory, "windows")
    data_info = join(inversion_directory, "data_info")
    stations = join(inversion_directory, "stations")
    raw_sync = join(inversion_directory, "syncs", "raw_sync")
    raw_cmtfiles = join(inversion_directory, "cmts", "raw_cmtsolutions")
    # directories for each iteration
    iter_green1_cmt = join(inversion_directory, "cmts",
                           f"iter{iter_number}_green1_cmtsolutions")
    iter_green1_sync = join(inversion_directory, "syncs",
                            f"iter{iter_number}_green1_sync")
    iter_conv1_processed = join(
        inversion_directory, "syncs", f"iter{iter_number}_conv1_processed_sync")
    iter_misfit_windows = join(
        inversion_directory, "misfit_windows", f"iter{iter_number}_misfit_windows")
    iter_adjoint_source = join(
        inversion_directory, "adjoint_sources", f"iter{iter_number}_adjoint_sources")
    iter_src_frechets = join(
        inversion_directory, "src_frechets", f"iter{iter_number}_src_frechets")
    iter_green2_cmt = join(inversion_directory, "cmts",
                           f"iter{iter_number}_green2_cmtsolutions")
    iter_green2_sync = join(inversion_directory, "syncs",
                            f"iter{iter_number}_green2_sync")
    iter_next_cmt = join(inversion_directory, "cmts",
                         f"iter{iter_number}_next_cmtsolutions")

    # make directories
    if (int(iter_number) == 1):
        sh.mkdir("-p", specfem_base)
        sh.mkdir("-p", specfem_stations_adjoint)
        sh.mkdir("-p", join(inversion_directory, "cmts"))
        sh.mkdir("-p", join(inversion_directory, "syncs"))
        sh.mv(stations_path, stations)
        sh.mv(cmtfiles_directory, raw_cmtfiles)
        sh.cp("-r", raw_cmtfiles, specfem_cmtfiles)
        sh.mv(ref_directory, specfem_ref)
        sh.mkdir("-p", specfem_output)
        sh.mkdir("-p", specfem_database)
        sh.mv(data_directory, data_processed)
        sh.mv(windows_directory, windows)
        sh.mv(data_info_directory, data_info)
        sh.mv(raw_sync_directory, raw_sync)
    elif (int(iter_number) > 1):
        sh.rm("-rf", specfem_cmtfiles)
        iter_last_cmt = join(inversion_directory, "cmts",
                             f"iter{int(iter_number)-1}_next_cmtsolutions")
        sh.cp("-r", iter_last_cmt, specfem_cmtfiles)

    sh.mkdir("-p", iter_green1_cmt)
    sh.mkdir("-p", iter_green1_sync)
    sh.mkdir("-p", iter_conv1_processed)
    sh.mkdir("-p", iter_misfit_windows)
    sh.mkdir("-p", iter_adjoint_source)
    sh.mkdir("-p", iter_src_frechets)
    sh.mkdir("-p", iter_green2_cmt)
    sh.mkdir("-p", iter_green2_sync)
    sh.mkdir("-p", iter_next_cmt)
    return (specfem_base, specfem_stations_adjoint, specfem_cmtfiles, specfem_ref, specfem_output, specfem_database, data_processed, windows, data_info, stations, raw_sync, iter_green1_cmt, iter_green1_sync,
            iter_conv1_processed, iter_misfit_windows, iter_adjoint_source, iter_src_frechets, iter_green2_cmt,
            iter_green2_sync, iter_next_cmt)


@click.command()
@click.option('--iter_number', required=True, type=int, help="the iteration number, the first one is 1")
@click.option('--py', required=True, type=str, help="the python path")
@click.option('--n_total', required=True, type=int, help="total number of events")
@click.option('--n_each', required=True, type=int, help="number of events for each iteration in specfem simulation")
@click.option('--n_iter', required=True, type=int, help="number of iterations in specfem simulation")
@click.option('--nproc', required=True, type=int, help="the number of processes for each specfem simulation")
@click.option('--n_node', required=True, type=int, help="the total number of nodes to use")
@click.option('--ntasks', required=True, type=int, help="the total number of mpi processes to use")
@click.option('--partition', required=True, type=str, help="the partion to use")
@click.option('--simulation_time_step1', required=True, type=str, help="the simulation time for step 1")
@click.option('--account', required=True, type=str, help="the account to use")
@click.option('--n_node_line_search', required=True, type=int, help="the total number of nodes to use in step 2")
@click.option('--ntasks_line_search', required=True, type=int, help="the total number of processes to use in step 2")
@click.option('--partition_line_search', required=True, type=str, help="the partion to use in step 2")
@click.option('--simulation_time_step2', required=True, type=str, help="the simulation time for step 2")
@click.option('--inversion_directory', required=True, type=str, help="the inversion directory")
@click.option('--cmtfiles_directory', required=True, type=str, help="the cmtfiles directory")
@click.option('--ref_directory', required=True, type=str, help="the specfem reference directory")
@click.option('--data_directory', required=True, type=str, help="the processed data directory")
@click.option('--windows_directory', required=True, type=str, help="the windows directory")
@click.option('--data_info_directory', required=True, type=str, help="the datainfo directory")
@click.option('--stations_path', required=True, type=str, help="the stations path")
@click.option('--raw_sync_directory', required=True, type=str, help="the raw sync directory (not changing tau)")
@click.option('--waveform_length', required=True, type=int, help="the waveform length")
@click.option('--taper_tmin_tmaxs', required=True, type=str, help="the taper time bands: minp1,maxp1/minp2,maxp2/...")
@click.option('--periods', required=True, type=str, help="min periods in filtering: minp1,maxp1/minp2,maxp2/...")
@click.option('--sampling_rate', required=True, type=int, help="the sampling rate to use")
@click.option('--alpha_range', required=True, type=str, help="the line search range for alpha")
@click.option('--t0_range', required=True, type=str, help="the line search range for t0")
@click.option('--tau_range', required=True, type=str, help="the line search range for tau, use 'fixed' if no change")
def source_inversion_single_step(iter_number, py, n_total, n_each, n_iter, nproc, n_node, ntasks, partition, simulation_time_step1, account,
                                 n_node_line_search, ntasks_line_search, partition_line_search, simulation_time_step2,
                                 inversion_directory, cmtfiles_directory, ref_directory, data_directory, windows_directory, data_info_directory, stations_path,
                                 raw_sync_directory,
                                 waveform_length, taper_tmin_tmaxs, periods, sampling_rate,
                                 alpha_range, t0_range, tau_range):
    """
    perform the source inversion for one step.
    """
    # firstly we have to make a running directory.
    (specfem_base, specfem_stations_adjoint, specfem_cmtfiles, specfem_ref, specfem_output, specfem_database, data_processed, windows, data_info, stations, raw_sync, iter_green1_cmt, iter_green1_sync,
     iter_conv1_processed, iter_misfit_windows, iter_adjoint_source, iter_src_frechets, iter_green2_cmt,
     iter_green2_sync, iter_next_cmt) = make_source_inversion_directory(
        iter_number, inversion_directory, cmtfiles_directory, ref_directory, data_directory, windows_directory, data_info_directory, stations_path,
        raw_sync_directory)

    # * step 1: generaet two green function asdf files.
    # init
    current_directory = str(sh.pwd())[:-1]
    result = "date; "
    result += "module load boost/1.68; "
    result += "module load phdf5/1.8.16; "
    # later we generate the cmt solution files for the green function
    result += generate_green_cmtsolutions(py,
                                          specfem_cmtfiles, iter_green1_cmt)
    # make simulation dirs based on the green function
    if(int(iter_number) == 1):
        init_structure(specfem_base, specfem_cmtfiles, specfem_ref,
                       specfem_output, specfem_database)
    # here we generate the directory using specfem_cmtfiles, but have to change to iter_green1_cmt
    result += cp_cmtsolution2structure(py,
                                       iter_green1_cmt, specfem_base)
    # change simulation type to forward
    result += change_simulation_type(py, specfem_base, "forward")
    # do forward simulation for green1
    result += forward_task(base=specfem_base, N_total=n_total,
                           N_each=n_each, N_iter=n_iter, nproc=nproc, run_mesh=True)
    result += f"cd {current_directory}; \n"
    # collect the sync to green1
    result += collect_sync_files(py, specfem_output, iter_green1_sync)
    # convert green1 to conv1
    result += convert_green2conv_processed(py, n_total, specfem_cmtfiles, iter_green1_sync,
                                           iter_conv1_processed, waveform_length, taper_tmin_tmaxs, periods, sampling_rate)
    # calculate the misfit windows
    body_periods, surface_periods = periods.split("/")
    body_min_period, body_max_period = body_periods.split(",")
    surface_min_period, surface_max_period = surface_periods.split(",")
    min_periods = f"{body_min_period},{surface_min_period}"
    max_periods = f"{body_max_period},{surface_max_period}"
    result += calculate_misfit_windows(py, n_total, windows, iter_misfit_windows,
                                       min_periods, max_periods, data_processed, iter_conv1_processed, data_info)
    # generate stations_adjoint directory and copy it to the data directory
    result += calculate_stations_adjoint(py, stations,
                                         iter_misfit_windows, specfem_stations_adjoint)
    result += cp_stations_adjoint2structure(py,
                                            specfem_stations_adjoint, specfem_base)
    # calculate the adjoint source
    result += calculate_adjoint_source_raw(py, n_total, iter_misfit_windows, stations, raw_sync, iter_conv1_processed,
                                           data_processed, data_info, iter_adjoint_source, body_periods, surface_periods)
    # move the adjoint source back to the SEM folder
    result += ln_adjoint_source_to_structure(py,
                                             iter_adjoint_source, specfem_base)
    # change simulation type to the source inversion (type 2)
    result += change_simulation_type(py, specfem_base, "source")
    # do the adjoint simulation
    result += forward_task(base=specfem_base, N_total=n_total,
                           N_each=n_each, N_iter=n_iter, nproc=nproc, run_mesh=False)
    result += f"cd {current_directory}; \n"
    # collect the src_frechet files
    result += collect_src_frechet_files(py, specfem_output, iter_src_frechets)
    # make the perturbed green cmtsolutions
    result += make_perturbed_cmtsolution(py, iter_src_frechets,
                                         iter_green1_cmt, iter_green2_cmt)
    # cp the perturbed cmtsolution to the simulation directory
    result += cp_cmtsolution2structure(py,
                                       iter_green2_cmt, specfem_base)
    # change simulation type to forward
    result += change_simulation_type(py, specfem_base, "forward")
    # do the forward simulation
    result += forward_task(base=specfem_base, N_total=n_total,
                           N_each=n_each, N_iter=n_iter, nproc=nproc, run_mesh=False)
    result += f"cd {current_directory}; \n"
    # collect the sync to green2
    result += collect_sync_files(py, specfem_output, iter_green2_sync)

    # * now we can submit the job
    cal_asdf_job_id = submit_job(
        f"iter{iter_number}_step1_source", result, n_node, ntasks, partition, simulation_time_step1, account, "stampede2")

    # * write the job script to perform the line search for this step
    # note if the first job finished before submitting the second job, there will be problem. (appearantly not the case)
    # init
    result = "date; "
    result += "module load boost/1.68; "
    result += "module load phdf5/1.8.16; "
    # do the source line search
    result += source_line_search(py, n_total, iter_green1_sync, iter_green2_sync, data_processed, windows, data_info, stations,
                                 iter_src_frechets, specfem_cmtfiles, iter_next_cmt,
                                 body_periods, surface_periods, taper_tmin_tmaxs,
                                 alpha_range, t0_range, tau_range)
    # * submit the line search job
    line_search_job_id = submit_job(
        f"iter{iter_number}_step2_source", result, n_node_line_search, ntasks_line_search, partition_line_search,
        simulation_time_step2, account, "stampede2", depends_on=[cal_asdf_job_id])


if __name__ == "__main__":
    source_inversion_single_step()
