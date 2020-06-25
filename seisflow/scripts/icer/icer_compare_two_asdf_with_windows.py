"""
icer_compare_two_asdf_with_windows.py: 
"""
from ...slurm.submit_job import submit_job
import click
import sys


def get_scripts(misfit_windows_dir, obs_dir, syn_dir, output_dir, info_dir, azimuth_width, waves_perpage, snr, cc, deltat, band, plot_surface, n_cores, py):
    result = ""
    result += "module purge;"
    result += "module load GCC/8.2.0-2.31.1;"
    result += "module load OpenMPI/3.1.3;"
    if(plot_surface):
        result += f"srun -n {n_cores} {py} -m seisflow.scripts.plot.mpi_compare_two_asdf_with_windows --misfit_windows_dir {misfit_windows_dir} --obs_dir {obs_dir} --syn_dir {syn_dir} --output_dir {output_dir} --info_dir {info_dir} --azimuth_width {azimuth_width} --waves_perpage {waves_perpage} --snr {snr} --cc {cc} --deltat {deltat} --band {band} --plot_surface ;"
    else:
        result += f"srun -n {n_cores} {py} -m seisflow.scripts.plot.mpi_compare_two_asdf_with_windows --misfit_windows_dir {misfit_windows_dir} --obs_dir {obs_dir} --syn_dir {syn_dir} --output_dir {output_dir} --info_dir {info_dir} --azimuth_width {azimuth_width} --waves_perpage {waves_perpage} --snr {snr} --cc {cc} --deltat {deltat} --band {band};"
    return result


@click.command()
@click.option('--misfit_windows_dir', required=True, type=str)
@click.option('--obs_dir', required=True, type=str)
@click.option('--syn_dir', required=True, type=str)
@click.option('--output_dir', required=True, type=str)
@click.option('--info_dir', required=True, type=str)
@click.option('--azimuth_width', required=True, type=int)
@click.option('--waves_perpage', required=True, type=int)
@click.option('--snr', required=True, type=float)
@click.option('--cc', required=True, type=float)
@click.option('--deltat', required=True, type=float)
@click.option('--band', required=True, type=str)
@click.option('--plot_surface', is_flag=True, default=False)
@click.option('--n_cores', required=True, type=int)
@click.option('--used_time', required=True, type=str)
def main(misfit_windows_dir, obs_dir, syn_dir, output_dir, info_dir, azimuth_width, waves_perpage, snr, cc, deltat, band, plot_surface, n_cores, used_time):
    py = sys.executable
    to_submit_script = get_scripts(misfit_windows_dir, obs_dir, syn_dir, output_dir, info_dir,
                                   azimuth_width, waves_perpage, snr, cc, deltat, band, plot_surface, n_cores, py)
    submit_job("traces_figure", to_submit_script, None,
               n_cores, None, used_time, None, "icer_flexiable")


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
