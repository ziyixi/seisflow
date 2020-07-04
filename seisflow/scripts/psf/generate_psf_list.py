"""
generate_psf_list.py: generate the psf list used by xsem_make_gauss_point_from_list.
"""
import click
import numpy as np


def lld2xyz(lon, lat, dep):
    R_earth = 6371.0
    r = (R_earth - dep) / R_earth
    theta = 90 - lat
    phi = lon

    z = r * np.cos(np.deg2rad(theta))
    h = r * np.sin(np.deg2rad(theta))
    x = h * np.cos(np.deg2rad(phi))
    y = h * np.sin(np.deg2rad(phi))

    return x, y, z


@click.command()
@click.option('--lons', required=True, type=str, help="lon1,lon2,steplon")
@click.option('--lats', required=True, type=str, help="lat1,lat2,steplat")
@click.option('--deps', required=True, type=str, help="dep1,dep2,stepdep")
@click.option('--gauss_width_km', required=True, type=float, help="the gaussian width in km")
@click.option('--perturbation', required=True, type=float, help="abs max perturbation")
@click.option('--output_path', required=True, type=str, help="the output path")
def main(lons, lats, deps, gauss_width_km, perturbation, output_path):
    lon1, lon2, steplon = map(float, lons.split(","))
    lat1, lat2, steplat = map(float, lats.split(","))
    dep1, dep2, stepdep = map(float, deps.split(","))
    used_lons = np.arange(lon1, lon2 + steplon / 2, steplon)
    used_lats = np.arange(lat1, lat2 + steplat / 2, steplat)
    used_deps = np.arange(dep1, dep2 + stepdep / 2, stepdep)
    # * firstly we generate the pos/neg perturbation value list (3D)
    # 2D first
    start_lat = -1
    plot_color = None
    colorlist_2d = np.zeros((len(used_lons), len(used_lats)))
    for ilat, _ in enumerate(used_lats):
        start_lat = -start_lat
        plot_color = -start_lat
        for ilon, _ in enumerate(used_lons):
            colorlist_2d[ilon, ilat] = plot_color
    # 3D
    colorlist_3d = np.zeros((len(used_lons), len(used_lats), len(used_deps)))
    for idep, _ in enumerate(used_deps):
        if (idep % 2 == 0):
            colorlist_3d[:, :, idep] = colorlist_2d
        else:
            colorlist_3d[:, :, idep] = -colorlist_2d

    # * now we can output the list
    with open(output_path, "w") as f:
        for ilon, each_lon in enumerate(used_lons):
            for ilat, each_lat in enumerate(used_lats):
                for idep, each_dep in enumerate(used_deps):
                    x, y, z = lld2xyz(each_lon, each_lat, each_dep)
                    model_value = colorlist_3d[ilon, ilat, idep] * perturbation
                    f.write(f"{x} {y} {z} {gauss_width_km} {model_value} \n")


if __name__ == "__main__":
    main()
