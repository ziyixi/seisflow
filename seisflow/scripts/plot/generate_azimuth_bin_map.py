"""
generate_azimuth_bin_map.py: give the gcmtid, plot the corresponding azimuth bins.
"""
from os.path import basename, dirname, join

import click
import numpy as np
import pygmt
import sh
import tqdm
import xarray as xr
from PyPDF2 import PdfFileMerger

from .generate_model_h_slices import plot_single_figure

# countries/region in Asia, except CN,TW,HK,MO,IN
countries = ["AF", "AM", "AZ", "BH", "BD", "BN", "KH", "CX", "CC", "IO", "GE", "ID", "IR", "IQ", "IL", "JP", "JO", "KZ", "KW", "KG", "LA", "LB",
             "MY", "MV", "MN", "NP", "KP", "OM", "PK", "PS", "PH", "QA", "SA", "SG", "KR", "LK", "SY", "TJ", "TH", "TR", "TM", "AE", "UZ", "VN", "YE"]
countries = ",".join(countries)


@click.command()
@click.option('--model_file', required=True, type=str, help="the model netcdf file path")
@click.option('--region', required=True, type=str, help="the region to plot, lon1/lat1/lon2/lat2")
@click.option('--npts', required=False, default="1000/1000", type=str, help="the npts as lonnpts/latnpts and also hnpts/vnpts")
@click.option('--vmin', required=False, default="-0.08", type=float, help="colorbar vmin")
@click.option('--vmax', required=False, default="0.08", type=float, help="colorbar vmax")
@click.option('--output_path', required=True, type=str, help="the output pdf path")
@click.option('--colorbar', required=False, default="seis", type=str, help="the colorbar")
@click.option('--paths', required=True, type=str, help="the paths file, each row: lon1 lat1 lon2 lat2 lon/lat")
def main(model_file, region, npts, vmin, vmax, output_path, colorbar, paths):
    data = xr.open_dataset(model_file)
    # * make a hidden dir in the output_path
    lonnpts, latnpts = map(int, npts.split("/"))
    lon1, lat1, lon2, lat2 = map(float, region.split("/"))
    # * we will only generate one figure
