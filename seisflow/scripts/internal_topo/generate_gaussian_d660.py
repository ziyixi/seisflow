"""
generate_gaussian_d660.py: generate DATA/s362ani/external_topo_410_660.txt in Specfem while using external internal topo.
"""
import click
import numpy as np


def generate_650_topo(latmin, latmax, dep0, lon0, lat0, R, output_path, lonmin, lonmax):
    if(lonmin == None and lonmax == None):
        lon1 = np.arange(0, 360.2, 0.2)
    else:
        lon1 = np.arange(lonmin, lonmax+0.2, 0.2)
    lat1 = np.arange(latmin, latmax+0.2, 0.2)
    lon2, lat2 = np.meshgrid(lon1, lat1, indexing="ij")
    one_sigma2 = R**2
    d650 = np.zeros_like(lon2)
    for ilat in range(len(lat1)):
        for ilon in range(len(lon1)):
            lat = lat2[ilon, ilat]
            lon = lon2[ilon, ilat]
            if((lon > 180) and (lonmin == None and lonmax == None)):
                lon = lon-360
            r = (lat-lat0)**2+(lon-lon0)**2
            d650[ilon, ilat] = 0 + dep0 * np.exp(-0.5 * r / one_sigma2)
    if(output_path != None):
        with open(output_path, "w") as f:
            f.write(f"{d650.shape[0]} {d650.shape[1]} \n")
            for ilat in range(len(lat1)):
                for ilon in range(len(lon1)):
                    f.write(
                        f"{lon2[ilon,ilat]:.2f} {lat2[ilon,ilat]:.2f} 0.0000 {d650[ilon,ilat]:.4f} \n")
    return d650


@click.command()
@click.option('--lat_range', required=True, type=str, help="latmin,latmax should be larger than the simulation region")
@click.option('--dep0', required=True, type=float, help="the depression(+) or uplifting(-) depth")
@click.option('--lon0', required=True, type=float, help="the longitude for the depression/uplifting")
@click.option('--lat0', required=True, type=float, help="the latitude for the depression/uplifting")
@click.option('--radius', required=True, type=float, help="the radius for the depression/uplifting")
@click.option('--output_path', required=True, type=str, help="the output path for the depression/uplifting file")
def main(lat_range, dep0, lon0, lat0, radius, output_path):
    latmin, latmax = map(float, lat_range.split(","))
    generate_650_topo(latmin, latmax, dep0, lon0, lat0,
                      radius, output_path, None, None)


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
