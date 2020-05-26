"""
generate_model_v_slices.py: generate the vertical slices for the model.
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
@click.option('--smooth_index', required=False, default="41,65", type=str, help="the depth index to smooth")
@click.option('--parameters', required=False, default="vsv,vsh,vpv,vph", type=str, help="the parameters to use")
@click.option('--vmin', required=False, default="-0.08", type=float, help="colorbar vmin")
@click.option('--vmax', required=False, default="0.08", type=float, help="colorbar vmax")
@click.option('--output_path', required=True, type=str, help="the output pdf path")
@click.option('--colorbar', required=False, default="seis", type=str, help="the colorbar")
@click.option('--paths', required=True, type=str, help="the paths file, each row: lon1 lat1 lon2 lat2 lon/lat")
def main(model_file, region, npts, smooth_index, parameters, vmin, vmax, output_path, colorbar, paths):
    data = xr.open_dataset(model_file)
    # * make a hidden dir in the output_path
    temp_directory = join(dirname(output_path), "."+basename(output_path))
    sh.mkdir("-p", temp_directory)
    # * some configurations
    parameters = parameters.split(",")
    lonnpts, latnpts = map(int, npts.split("/"))
    smooth_index = list(map(int, smooth_index.split(",")))
    lon1, lat1, lon2, lat2 = map(float, region.split("/"))
    pdfs = []

    # * firstly, we should generate a figure showing the paths we want to plot.
    hlat = np.linspace(lat1, lat2, latnpts)
    hlat = xr.DataArray(hlat, dims='hlat', coords={'hlat': hlat})
    hlon = np.linspace(lon1, lon2, lonnpts)
    hlon = xr.DataArray(hlon, dims='hlon', coords={'hlon': hlon})
    # * now we load the paths
    # * lon1 lat1 lon2 lat2 dep1 dep2 lon/lat
    paths = np.loadtxt(paths, dtype=np.str)
    plot_paths = []
    for row in paths:
        plot_paths.append([float(row[0]), float(row[1]),
                           float(row[2]), float(row[3])])
    # tqdm
    pbar = tqdm.tqdm(total=len(parameters)*len(plot_paths)+1)
    # use vsv to represent the paths
    to_interp_data = data["vsv"].copy()
    to_interp_data.data[to_interp_data.data > 9e6] = np.nan
    pdf_path_h_slice = plot_single_figure(to_interp_data, 100, hlat, hlon, colorbar,
                                          vmin, vmax, lon1, lon2, lat1, lat2, "vsv", temp_directory, plot_paths=plot_paths)
    pdfs.append(pdf_path_h_slice)
    pbar.update(1)
    # * now we plot vertical slices accordingly
    for line_index, row in enumerate(paths):
        # the bellow loop will not be so expensive.
        for each_parameter in parameters:
            to_interp_data = data[each_parameter].copy()
            for each_smooth_index in smooth_index:
                to_interp_data[:, :, each_smooth_index].data[:] = (
                    to_interp_data[:, :, each_smooth_index - 1].data + to_interp_data[:, :, each_smooth_index + 1].data) / 2
            to_interp_data.data[to_interp_data.data > 9e6] = np.nan
            pdf_path = plot_single_vertical_figure(
                line_index, row, to_interp_data, npts, each_parameter, colorbar, vmin, vmax, temp_directory)
            pdfs.append(pdf_path)
            pbar.update(1)

    merger = PdfFileMerger()
    for pdf in pdfs:
        merger.append(pdf)
    merger.write(output_path)
    merger.close()
    pbar.close()


def plot_single_vertical_figure(line_index, each_path, to_interp_data, npts, each_parameter, colorbar, vmin, vmax, temp_directory):
    # * firstly, we generate h and v, two new axis.
    hnpts, vnpts = map(float, npts.split("/"))
    lon1, lat1, lon2, lat2, dep1, dep2, hlabel = each_path
    lon1, lat1, lon2, lat2, dep1, dep2 = map(
        float, [lon1, lat1, lon2, lat2, dep1, dep2])
    newlat = np.linspace(lat1, lat2, hnpts)
    newlon = np.linspace(lon1, lon2, hnpts)
    newdep = np.linspace(dep1, dep2, vnpts)
    if (hlabel == "lon"):
        hrep = newlon
    else:
        hrep = newlat
    newlat_xr = xr.DataArray(newlat, dims='h', coords={'h': hrep})
    newlon_xr = xr.DataArray(newlon, dims='h', coords={'h': hrep})
    newdep_xr = xr.DataArray(newdep, dims='v', coords={'v': newdep})
    # * generate the plotting data
    plot_data = to_interp_data.interp(
        depth=newdep_xr, latitude=newlat_xr, longitude=newlon_xr)
    plot_data.coords["v"] = 6371 - plot_data.coords["v"]
    plot_data = plot_data.T
    # * now we can plot the figure
    fig = pygmt.Figure()
    if (colorbar == "default"):
        pygmt.makecpt(cmap="seisflow/data/dvs_6p.cpt", series=f"{vmin}/{vmax}/0.01",
                      continuous=True, D="o")
    elif(colorbar[-2:] == "_r"):
        pygmt.makecpt(cmap=colorbar[:-2], series=f"{vmin}/{vmax}/0.01",
                      continuous=True, D="o", reverse=True)
    else:
        pygmt.makecpt(cmap=colorbar, series=f"{vmin}/{vmax}/0.01",
                      continuous=True, D="o")
    thetitle = f"Line{line_index},{each_parameter}"
    if(hlabel == "lon"):
        fig.basemap(projection=f"Pa20c/{(lon1+lon2)/2}z",
                    region=f"{lon1}/{lon2}/{6371-dep2}/6371", frame=["afg", f"+t{thetitle}", "+lLontitude"])
    else:
        fig.basemap(projection=f"Pa20c/{(lat1+lat2)/2}z",
                    region=f"{lat1}/{lat2}/{6371-dep2}/6371", frame=["afg", f"+t{thetitle}", "+lLatitude"])
    fig.grdimage(plot_data, cmap=True)
    fig.colorbar(
        # justified inside map frame (j) at Top Center (TC)
        position="JBC+w14c/1.5c+h+e",
        box=False,
        frame=[f"+LdlnV{each_parameter[1:]}", "xaf"],
        scale=100,)
    pdf_path = join(
        temp_directory, f"{each_parameter}_{lon1},{lat1}_{lon2},{lat2}.pdf")
    fig.savefig(
        pdf_path)
    return pdf_path


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
