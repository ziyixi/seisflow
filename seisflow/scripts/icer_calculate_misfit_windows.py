"""
icer_calculate_misfit_windows.py: calculate misfit windows massively on ICER.
"""
from os import listdir

from ..slurm import submit_job


def get_script(windows_directory, output_directory, time_length, station_fname, min_periods, max_periods, data_asdf_directory, sync_asdf_directory, data_info_directory,
               n_node, n_cores_each_node, py, used_time):
    """
    get_script: get slurm script.
    """
    result = ""
    result += "module purge;"
    result += "module load GCC/8.2.0-2.31.1;"
    result += "module load OpenMPI/3.1.3;"
    # run mpi_extract_data_info.py
    if (len(listdir(data_info_directory)) == 0):
        # run only when data_info_directory is empty
        result += f"srun -n {n_node*n_cores_each_node} {py} -m seisflow.scripts.mpi_extract_data_info --asdf_directory {data_asdf_directory} --station_fname {station_fname} --output_dir {data_info_directory} ;"
    else:
        pass
    # get windows
    result += f"srun -n {n_node*n_cores_each_node} {py} -m seisflow.scripts.mpi_tao_2018_ggg_windows --data_info_directory {data_info_directory} --time_length {time_length} --output_dir {windows_directory} ;"
    # run mpi_calculate_misfit_windows.py
    result += f"srun -n {n_node*n_cores_each_node} {py} -m seisflow.scripts.mpi_calculate_misfit_windows --windows_directory {windows_directory} --output_directory {output_directory} --min_periods {min_periods} --max_periods {max_periods} --data_asdf_directory {data_asdf_directory} --sync_asdf_directory {sync_asdf_directory} --data_info_directory {data_info_directory} ;"
    return result


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--windows_directory', required=True, type=str, help="the directory to store windows")
    @click.option('--output_directory', required=True, type=str, help="the directory to output the misfit windows")
    @click.option('--data_asdf_directory', required=True, type=str, help="the data asdf directory")
    @click.option('--sync_asdf_directory', required=True, type=str, help="the sync asdf directory")
    @click.option('--data_info_directory', required=True, type=str, help="the data info directory to store data info")
    @click.option('--time_length', required=True, type=int, help="the time length used")
    @click.option('--station_fname', required=True, type=str, help="the stations path in Specfem format")
    @click.option('--min_periods', required=True, type=str, help="the min periods in filtering: min_body,min_surface")
    @click.option('--max_periods', required=True, type=str, help="the max periods in filtering: max_body,max_surface")
    @click.option('--n_node', required=True, type=int, help="the number of nodes to use")
    @click.option('--n_cores_each_node', required=True, type=int, help="the number of cores to use for each node")
    @click.option('--py', required=True, type=str, help="the python path")
    @click.option('--used_time', required=True, type=str, help="the time to use")
    def main(windows_directory, output_directory, time_length, station_fname, min_periods, max_periods, data_asdf_directory, sync_asdf_directory, data_info_directory,
             n_node, n_cores_each_node, py, used_time):
        to_submit_script = get_script(windows_directory, output_directory, time_length, station_fname, min_periods, max_periods, data_asdf_directory, sync_asdf_directory, data_info_directory,
                                      n_node, n_cores_each_node, py, used_time)
        submit_job("calculate_misfit", to_submit_script, n_node, n_node *
                   n_cores_each_node, None, used_time, None, "icer")
    main()
