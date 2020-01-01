"""
process_data.py: the wrapper to process data for one asdf file.
"""
from ..asdf.process_data import process_single_event

if __name__ == "__main__":
    import click
    from .. import get_data

    @click.command()
    @click.option('--min_periods', required=True, type=str, help="min periods in seconds, eg: 10,40")
    @click.option('--max_periods', required=True, type=str, help="max periods in seconds, eg: 120,120")
    @click.option('--taper_tmin_tmax', required=True, type=str, help="frequency taper tmin(f3) and tmax(f2), eg: 1,400")
    @click.option('--asdf_filename', required=True, type=str, help="asdf raw data file name")
    @click.option('--waveform_length', required=True, type=float, help="waveform length to cut (from event start time)")
    @click.option('--sampling_rate', required=True, type=int, help="sampling rate in HZ")
    @click.option('--output_directory', required=True, type=str, help="output directory")
    @click.option('--correct_cea/--no-correct_cea', default=False, help="if handling cea dataset")
    def main(min_periods, max_periods, taper_tmin_tmax, asdf_filename, waveform_length, sampling_rate, output_directory, logfile, correct_cea):
        min_periods = [float(item) for item in min_periods.split(",")]
        max_periods = [float(item) for item in max_periods.split(",")]
        cea_correction_file=get_data("cmpaz_segment.txt")
        process_single_event(min_periods, max_periods, taper_tmin_tmax, asdf_filename,
                             waveform_length, sampling_rate, output_directory, correct_cea, cea_correction_file)

    main()