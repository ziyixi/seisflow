"""
generate_model_h_slices.py: use pygmt to generate the model horizantal slices into a pdf file.
"""
from os.path import basename, dirname, join

import click
import numpy as np
import pygmt
import sh
import tqdm
import xarray as xr
from PyPDF2 import PdfFileMerger

# countries/region in Asia, except CN,TW,HK,MO,IN
countries = ["AF", "AM", "AZ", "BH", "BD", "BN", "KH", "CX", "CC", "IO", "GE", "ID", "IR", "IQ", "IL", "JP", "JO", "KZ", "KW", "KG", "LA", "LB",
             "MY", "MV", "MN", "NP", "KP", "OM", "PK", "PS", "PH", "QA", "SA", "SG", "KR", "LK", "SY", "TJ", "TH", "TR", "TM", "AE", "UZ", "VN", "YE"]
countries = ",".join(countries)


@click.command()
@click.option('--model_file', required=True, type=str, help="the model netcdf file path")
@click.option('--region', required=True, type=str, help="the region to plot, lon1/lat1/lon2/lat2")
@click.option('--npts', required=False, default="1000/1000", type=str, help="the npts as lonnpts/latnpts")
@click.option('--smooth_index', required=False, default="41,65", type=str, help="the depth index to smooth")
@click.option('--parameters', required=False, default="vsv,vsh,vpv,vph", type=str, help="the parameters to use")
@click.option('--depths', required=False, default="100,200,300,400,500", type=str, help="the depths to plot")
@click.option('--vmin', required=False, default="-0.08", type=float, help="colorbar vmin")
@click.option('--vmax', required=False, default="0.08", type=float, help="colorbar vmax")
@click.option('--output_path', required=True, type=str, help="the output pdf path")
@click.option('--colorbar', required=False, default="seis", type=str, help="the colorbar")
def main(model_file, region, npts, smooth_index, parameters, depths, vmin, vmax, output_path, colorbar):
    data = xr.open_dataset(model_file)
    # * make a hidden dir in the output_path
    temp_directory = join(dirname(output_path), "."+basename(output_path))
    sh.mkdir("-p", temp_directory)
    # * loop for each depth and each parameter
    parameters = parameters.split(",")
    depths = list(map(float, depths.split(",")))
    lonnpts, latnpts = map(int, npts.split("/"))
    lon1, lat1, lon2, lat2 = map(float, region.split("/"))
    smooth_index = list(map(int, smooth_index.split(",")))
    pdfs = []
    pbar = tqdm.tqdm(total=len(parameters)*len(depths))
    for each_parameter in parameters:
        hlat = np.linspace(lat1, lat2, latnpts)
        hlat = xr.DataArray(hlat, dims='hlat', coords={'hlat': hlat})
        hlon = np.linspace(lon1, lon2, lonnpts)
        hlon = xr.DataArray(hlon, dims='hlon', coords={'hlon': hlon})
        to_interp_data = data[each_parameter].copy()
        for each_smooth_index in smooth_index:
            to_interp_data[:, :, each_smooth_index].data[:] = (
                to_interp_data[:, :, each_smooth_index - 1].data + to_interp_data[:, :, each_smooth_index + 1].data) / 2
        to_interp_data.data[to_interp_data.data > 9e6] = np.nan
        # some parameters
        for each_depth in depths:
            pdf_path = plot_single_figure(to_interp_data, each_depth, hlat, hlon, colorbar,
                                          vmin, vmax, lon1, lon2, lat1, lat2, each_parameter, temp_directory)
            pdfs.append(pdf_path)
            pbar.update(1)
    pbar.close()
    merger = PdfFileMerger()
    for pdf in pdfs:
        merger.append(pdf)
    merger.write(output_path)
    merger.close()


def plot_single_figure(to_interp_data, each_depth, hlat, hlon, colorbar, vmin, vmax, lon1, lon2, lat1, lat2, each_parameter, temp_directory, plot_paths=None):
    plot_data = to_interp_data.interp(
        depth=each_depth, latitude=hlat, longitude=hlon)
    plot_data = plot_data.T
    # * plot for each depth and parameter
    fig = pygmt.Figure()
    if(colorbar[-2:] == "_r"):
        pygmt.makecpt(cmap=colorbar[:-2], series=f"{vmin}/{vmax}/0.01",
                      continuous=True, D="o", reverse=True)
    else:
        pygmt.makecpt(cmap=colorbar, series=f"{vmin}/{vmax}/0.01",
                      continuous=True, D="o")
    mean_lat = (lat1 + lat2) / 2
    mean_lon = (lon1 + lon2) / 2
    lat_dev = (max(lat1, lat2) - min(lat1, lat2)) / 10
    lon_dev = (max(lon1, lon2) - min(lon1, lon2)) / 10
    fig.basemap(projection=f"B{mean_lon}/{mean_lat}/{mean_lat-lat_dev}/{mean_lat+lat_dev}/40c",
                region=[lon1, lon2, lat1, lat2], frame=["afg"])
    fig.grdimage(plot_data, cmap=True)
    fig.coast(shorelines="1/0.5p,black",
              E=f"{countries}+p1/0.2p,black,dashed")
    # * plot the correct Chinese boundary
    # hard code the path, so the script need to be called at the root path
    data_path = "seisflow/data/CN-border-L1.dat"
    fig.plot(data=data_path,
             pen="1/0.2p,black,dashed")
    fig.colorbar(
        # justified inside map frame (j) at Top Center (TC)
        position="JBC+w30c/1.5c+h+e",
        box=False,
        frame=[f"+LdlnV{each_parameter[1:]}(%)", "xaf"],
        scale=100,)
    fig.text(x=[lon1 + lon_dev, lon1 + lon_dev * 5 / 4], y=[lat2 - lat_dev / 2, lat2-lat_dev], text=[
        each_parameter, f"{each_depth}km"], font="30p,4,white+jMC")
    # * sometimes we want to plot some lines in map
    if (plot_paths != None):
        for index, item in enumerate(plot_paths):
            fig.plot(x=[item[0], item[2]], y=[
                     item[1], item[3]], pen="1.5p,black")
            fig.text(x=item[0], y=item[1],
                     text=f"{index}", font="15p,4,black+jMC")
    pdf_path = join(
        temp_directory, f"{each_parameter}_{each_depth}km.pdf")
    fig.savefig(
        pdf_path)
    return pdf_path


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
