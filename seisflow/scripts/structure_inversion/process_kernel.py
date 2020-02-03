"""
process_kernel.py: submit a job to process kernel. Firstly do the kernel summation along with preconditioning, later smooth the kernel. In the last generate a new model with LINE_SEARCH_PERTURBATION perturbation for doing line search.
"""
import sys
from glob import glob
from os.path import basename, join

import sh

from ...slurm.submit_job import submit_job
from ...utils.setting import LINE_SEARCH_PERTURBATION


def construct_structure(database_directory, ref_directory, sem_utils_directory, kernel_process_directory, input_model_directory):
    """
    construct_structure: construct the structure used in structure inversion.
    """
    # * construct kernel_process_directory based on ref_directory
    sh.mkdir("-p", kernel_process_directory)
    sh.ln("-s", join(ref_directory, "bin"),
          join(kernel_process_directory, "bin"))
    sh.ln("-s", join(ref_directory, "DATA"),
          join(kernel_process_directory, "DATA"))
    # * construct the sem_utils_bin directory for sem_utils
    sh.ln("-s", join(sem_utils_directory, "bin"),
          join(kernel_process_directory, "sem_utils_bin"))
    # * kernels_list.txt
    all_database_paths = sorted(glob(join(database_directory, "*")))
    all_gcmtids = [basename(item) for item in all_database_paths]
    kernels_list_path = join(kernel_process_directory, "kernels_list.txt")
    with open(kernels_list_path, "w") as file:
        for each_path in all_gcmtids:
            file.write(f"{each_path}\n")
    # * link database directory
    to_link_input_kernels_path = join(
        kernel_process_directory, "INPUT_KERNELS")
    try:
        sh.unlink(to_link_input_kernels_path)
    except sh.ErrorReturnCode_1:
        pass
    sh.ln("-s", database_directory, to_link_input_kernels_path)
    # * make OUTPUT_SUM directory
    to_make_output_sum_path = join(kernel_process_directory, "OUTPUT_SUM")
    sh.mkdir("-p", to_make_output_sum_path)
    # * prepare directories for model updating
    # INPUT_MODEL should be the gll directory
    to_link_input_model_path = join(kernel_process_directory, "INPUT_MODEL")
    try:
        sh.unlink(to_link_input_model_path)
    except sh.ErrorReturnCode_1:
        pass
    sh.ln("-s", input_model_directory, to_link_input_model_path)
    # we should make a INPUT_GRADIENT directory
    sh.mkdir("-p", join(kernel_process_directory, "INPUT_GRADIENT"))
    # the topo directory should be provided as one database path
    rep_database_path = all_database_paths[0]
    to_link_topo_path = join(kernel_process_directory, "topo")
    try:
        sh.unlink(to_link_topo_path)
    except sh.ErrorReturnCode_1:
        pass
    sh.ln("-s", rep_database_path, to_link_topo_path)
    # link topo to DATABASES_MPI
    sh.ln("-s", rep_database_path, join(kernel_process_directory, "DATABASES_MPI"))
    # make OUTPUT_MODEL
    sh.mkdir("-p", join(kernel_process_directory, "OUTPUT_MODEL"))
    # make smoothed directory
    sh.mkdir("-p", join(kernel_process_directory, "SMOOTHED_KERNEL"))


def do_preconditioned_summation(kernel_process_directory):
    """
    do_preconditioned_summation: get the script to do the preconditioned summation.
    """
    result = ""
    current_path = str(sh.pwd())[:-1]  # pylint: disable=not-callable
    result += f"cd {kernel_process_directory};"
    result += f"ibrun ./bin/xsum_preconditioned_kernels;"
    result += f"cd {current_path};\n"
    return result


def do_smoothing(kernel_process_directory, sigma_h, sigma_v, input_dir, output_dir, n_tasks):
    """
    do_smoothing: perform smoothing for the summed kernel. (use the workflow order in our lab)
    """
    # * the commented part is for using smoother in specfem, which could be very slow
    # result = ""
    # to_smooth_kernel_names = [
    #     "bulk_c_kernel", "bulk_betav_kernel", "bulk_betah_kernel", "eta_kernel"]
    # current_path = str(sh.pwd())[:-1]  # pylint: disable=not-callable
    # result += f"cd {kernel_process_directory};"
    # for each_name in to_smooth_kernel_names:
    #     result += f"ibrun ./bin/xsmooth_sem {sigma_h} {sigma_v} {each_name} {input_dir} {output_dir};"
    # result += f"cd {current_path};\n"
    # * if we use the smoother in sem_utils
    result = ""
    to_smooth_kernel_names = [
        "bulk_c_kernel", "bulk_betav_kernel", "bulk_betah_kernel", "eta_kernel"]
    to_smooth_kernel_names = ",".join(to_smooth_kernel_names)
    current_path = str(sh.pwd())[:-1]  # pylint: disable=not-callable
    result += f"cd {kernel_process_directory};"
    result += f"ibrun ./sem_utils_bin/xsem_smooth {n_tasks} {join(kernel_process_directory,'topo')} {input_dir} {to_smooth_kernel_names} {sigma_h} {sigma_v} {output_dir} _smooth;"
    result += f"cd {current_path};\n"
    return result


def ln_smoothed_kernel_to_input_dir(pyexec, smooth_dir, kernel_process_directory):
    """
    ln_smoothed_kernel_to_input_dir: remove the flag smooth and link the smoothed kernels to the directory INPUT_GRADIENT.
    """
    result = ""
    current_path = str(sh.pwd())[:-1]  # pylint: disable=not-callable
    result += f"cd {current_path};"
    result += f"ibrun -n 1 {pyexec} -m seisflow.scripts.structure_inversion.ln_smoothed_kernel_to_input_dir --smooth_dir {smooth_dir} --kernel_process_directory {kernel_process_directory};"
    result += f"cd {current_path};\n"
    return result


def iter1_generate_perturbed_kernel(kernel_process_directory, perturbed_value):
    """
    iter1_generate_perturbed_kernel: generate the perturbed kernel used for line search. iter1 will use only the steepest descent method.
    """
    result = ""
    current_path = str(sh.pwd())[:-1]  # pylint: disable=not-callable
    result += f"cd {kernel_process_directory};"
    # when we run add_model_globe_tiso, INPUT_MODEL(from previous gll directory), INPUT_GRADIENT(from ln smooth), topo(from database) have all been established.
    result += f"ibrun ./bin/xadd_model_tiso {perturbed_value};"
    # we should move the kernel files in OUTPUT_MODEL to perturbed_{perturbed_value}_for_line_search
    result += f"ibrun -n 1 mkdir -p perturbed_{perturbed_value}_for_line_search;"
    result += f"ibrun -n 1 mv OUTPUT_MODEL/* perturbed_{perturbed_value}_for_line_search/;"
    result += f"cd {current_path};\n"
    return result


def update_model_given_step_length(kernel_process_directory, perturbed_value):
    """
    update_model_given_step_length: update the model by optimized step length.
    """
    result = ""
    current_path = str(sh.pwd())[:-1]  # pylint: disable=not-callable
    result += f"cd {kernel_process_directory};"
    # when we run add_model_globe_tiso, INPUT_MODEL(from previous gll directory), INPUT_GRADIENT(from ln smooth), topo(from database) have all been established.
    result += f"ibrun ./add_model_globe_tiso {perturbed_value};"
    # we should move the kernel files in OUTPUT_MODEL to perturbed_{perturbed_value}_for_next_iteration
    result += f"ibrun -n 1 mkdir -p perturbed_{perturbed_value}_for_next_iteration;"
    result += f"ibrun -n 1 mv OUTPUT_MODEL/* perturbed_{perturbed_value}_for_next_iteration/;"
    result += f"cd {current_path};\n"
    return result


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--database_directory', required=True, type=str, help="the database directory")
    @click.option('--ref_directory', required=True, type=str, help="the reference specfem directory")
    @click.option('--sem_utils_directory', required=True, type=str, help="the reference sem_utils directory")
    @click.option('--kernel_process_directory', required=True, type=str, help="the directory to be created to process the kernels")
    @click.option('--input_model_directory', required=True, type=str, help="the input model directory")
    @click.option('--sigma_h', required=True, type=float, help="the value of sigma_h (km)")
    @click.option('--sigma_v', required=True, type=float, help="the value of sigma_v (km)")
    @click.option('--n_node', required=True, type=int, help="the number of nodes to use")
    @click.option('--n_tasks', required=True, type=int, help="the number of tasks to run")
    @click.option('--partition', required=True, type=str, help="the partion name, eg: skx-normal")
    @click.option('--time', required=True, type=str, help="the time used for processing the kernels")
    @click.option('--account', required=True, type=str, help="the account used in stampede2")
    def main(database_directory, ref_directory, sem_utils_directory, kernel_process_directory, input_model_directory, sigma_h, sigma_v, n_node, n_tasks,
             partition, time, account):
        """
        The main program to do the processing manually.
        """
        construct_structure(database_directory, ref_directory, sem_utils_directory,
                            kernel_process_directory, input_model_directory)
        # * prepare job script and submit the job in the later steps.
        result = ""
        result = "date; "
        result += "module load boost/1.68; "
        result += "module load phdf5/1.8.16;\n"
        # sum kernels, hessians and do precondition
        result += do_preconditioned_summation(kernel_process_directory)
        # do smoothing
        input_smooth_dir = join(kernel_process_directory, "OUTPUT_SUM")
        output_smooth_dir = join(kernel_process_directory, "SMOOTHED_KERNEL")
        result += do_smoothing(kernel_process_directory,
                               sigma_h, sigma_v, input_smooth_dir, output_smooth_dir, n_tasks)
        # ln smoothed kernels to the input directory
        pyexec = sys.executable
        result += ln_smoothed_kernel_to_input_dir(
            pyexec, output_smooth_dir, kernel_process_directory)
        # generate perturbed kernels with LINE_SEARCH_PERTURBATION step length for doing line search
        result += iter1_generate_perturbed_kernel(
            kernel_process_directory, LINE_SEARCH_PERTURBATION)
        # * now we can submit the job
        submit_job("process_kernel", result, n_node, n_tasks,
                   partition, time, account, "stampede2", depends_on=None)
    main()  # pylint: disable=no-value-for-parameter
