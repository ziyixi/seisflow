"""
icer_process_sync.py: massively process the sync on ICER.
"""
import sys
from glob import glob
from os.path import basename, join

from ...slurm.submit_job import submit_job


def get_files(base_dir):
    """
    get_files: get all the asdf files in base_dir.
    """
    return sorted(glob(join(base_dir, "*h5")))


def get_scripts(run_files, N_iters, N_node, N_cores_each_node, PY, min_periods, max_periods, waveform_length,
                sampling_rate, PROCESSED_DIR, taper_tmin_tmax, paz_directory):
    N_files = len(run_files)
    result = ""
    result += "module purge;"
    result += "module load GCC/8.2.0-2.31.1;"
    result += "module load OpenMPI/3.1.3;"
    # run iters
    for iiter in range(N_iters):
        result += f"echo 'start iteration {iiter}'; "
        for ieach in range(N_node):
            # run N_node files at the same iter
            offset = iiter*N_node+ieach
            if(offset >= N_files):
                continue
            filename = run_files[offset]
            # get paz file path
            filename_basename = basename(filename)
            gcmtid = filename_basename.split(".")[0].split(
                "_")[-1]  # work for raw_[gcmtid].h5 and [gcmtid].h5
            paz_path = join(paz_directory, gcmtid, "PZ")
            result += f"srun -n {N_cores_each_node} --exclusive {PY} -m seisflow.scripts.asdf.process_data_pz --min_periods {min_periods} --max_periods {max_periods} --taper_tmin_tmax {taper_tmin_tmax} --asdf_filename {filename} --waveform_length {waveform_length} --sampling_rate {sampling_rate} --output_directory {PROCESSED_DIR} --correct_cea --paz_path {paz_path} &"
        result += f"wait; "
        result += f"echo 'end iteration {iiter}'; "

    return result


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--data_directory', required=True, type=str, help="the raw sync directory")
    @click.option('--output_directory', required=True, type=str, help="the processed sync directory")
    @click.option('--n_iters', required=True, type=int, help="iterations to run")
    @click.option('--n_node', required=True, type=int, help="number of nodes to be used")
    @click.option('--n_cores_each_node', required=True, type=int, help="number of cores used for each node")
    @click.option('--used_time', required=True, type=str, help="the time to be used")
    @click.option('--periods', required=True, type=str, help="min periods in filtering: minp1,maxp1/minp2,maxp2/...")
    @click.option('--waveform_length', required=True, type=str, help="the length of the waveform to cut")
    @click.option('--sampling_rate', required=True, type=int, help="the sampling rate to use")
    @click.option('--taper_tmin_tmaxs', required=True, type=str, help="the taper time bands: minp1,maxp1/minp2,maxp2/...")
    @click.option('--paz_directory', required=True, type=str, help="paz directory")
    def main(data_directory, output_directory,
             n_iters, n_node, n_cores_each_node, used_time,
             periods, waveform_length, sampling_rate, taper_tmin_tmaxs, paz_directory):
        py = sys.executable
        run_files = get_files(data_directory)
        all_scripts = []
        for each_taper_tmin_tmax, each_period in zip(taper_tmin_tmaxs.split("/"), periods.split("/")):
            min_period, max_period = each_period.split(",")
            # note here we seprate different frequencies
            scripts_to_run = get_scripts(run_files, n_iters, n_node, n_cores_each_node, py, min_period, max_period, waveform_length,
                                         sampling_rate, output_directory, each_taper_tmin_tmax, paz_directory)
            all_scripts.append(scripts_to_run)
        to_submit_script = " ".join(all_scripts)
        submit_job("process_data", to_submit_script, n_node, n_node *
                   n_cores_each_node, None, used_time, None, "icer")

    main()  # pylint: disable=no-value-for-parameter
