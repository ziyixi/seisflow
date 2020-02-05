"""
xsede_process_sync.py: process the synthetics on stampede2.
"""

import sys
from glob import glob
from os.path import join

from ..slurm import submit_job


def get_files(base_dir):
    """
    get_files: get all the asdf files in base_dir.
    """
    return sorted(glob(join(base_dir, "*h5")))


def get_scripts(run_files, N_iters, n_event_each_iteration, N_cores_each_event, PY, min_periods, max_periods, waveform_length,
                sampling_rate, PROCESSED_DIR, taper_tmin_tmax):
    N_files = len(run_files)
    result = "date; "
    result += "module load boost/1.68;"
    result += "module load phdf5/1.8.16;"
    for iiter in range(N_iters):
        result += f"echo 'start iteration {iiter}'; "
        result += "date; "
        for ieach in range(n_event_each_iteration):
            # run N_node files at the same iter
            offset = iiter*n_event_each_iteration+ieach
            if(offset >= N_files):
                continue
            filename = run_files[offset]
            inc = ieach*N_cores_each_event
            result += f"ibrun -n {N_cores_each_event} -o {inc} {PY} -m seisflow.scripts.process_sync --min_periods {min_periods} --max_periods {max_periods} --asdf_filename {filename} --waveform_length {waveform_length} --sampling_rate {sampling_rate} --output_directory {PROCESSED_DIR} --taper_tmin_tmax {taper_tmin_tmax} &"
        result += f"wait; "
        result += f"echo 'end iteration {iiter}'; "
    return result


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--sync_directory', required=True, type=str, help="the raw sync directory")
    @click.option('--output_directory', required=True, type=str, help="the processed sync directory")
    @click.option('--n_iters', required=True, type=int, help="iterations to run")
    @click.option('--n_node', required=True, type=int, help="number of nodes to be used")
    @click.option('--n_event_each_iteration', required=True, type=int, help="number of events to process each iteration")
    @click.option('--n_cores_each_event', required=True, type=int, help="number of cores used for each node")
    @click.option('--used_time', required=True, type=str, help="the time to be used")
    @click.option('--periods', required=True, type=str, help="min periods in filtering: minp1,maxp1/minp2,maxp2/...")
    @click.option('--waveform_length', required=True, type=str, help="the length of the waveform to cut")
    @click.option('--sampling_rate', required=True, type=int, help="the sampling rate to use")
    @click.option('--taper_tmin_tmaxs', required=True, type=str, help="the taper time bands: minp1,maxp1/minp2,maxp2/...")
    @click.option('--partition', required=True, type=str, help="the partion to use, eg: skx-dev")
    @click.option('--account', required=True, type=str, help="the account in stampede2")
    def main(sync_directory, output_directory,
             n_iters, n_node, n_event_each_iteration, n_cores_each_event, used_time,
             periods, waveform_length, sampling_rate, taper_tmin_tmaxs, partition, account):
        run_files = get_files(sync_directory)
        all_scripts = []
        py = sys.executable
        for each_taper_tmin_tmax, each_period in zip(taper_tmin_tmaxs.split("/"), periods.split("/")):
            min_period, max_period = each_period.split(",")
            # note here we seprate different frequencies
            scripts_to_run = get_scripts(run_files, n_iters, n_event_each_iteration, n_cores_each_event, py, min_period, max_period, waveform_length,
                                         sampling_rate, output_directory, each_taper_tmin_tmax)
            all_scripts.append(scripts_to_run)
        to_submit_script = " ".join(all_scripts)
        # submit_job("process_sync", to_submit_script, n_node, n_event *
        #            n_cores_each_event, None, used_time, None, "icer")
        submit_job("process_sync", to_submit_script,
                   n_node, n_event_each_iteration * n_cores_each_event, partition, used_time, account, "stampede2")

    main()  # pylint: disable=no-value-for-parameter
