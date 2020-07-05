"""
generate_psf_list.py: generate the psf list used by xsem_make_gauss_point_from_list,
a wrap of https://github.com/taotaokai/sem_utils/blob/master/utils/point_spread_function/make_grid_xyz.py
"""
import click
import numpy as np
import pyproj


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
@click.option('--mesh_center_lat', required=True, type=float, help="the mesh center latitude")
@click.option('--mesh_center_lon', required=True, type=float, help="the mesh center longitude")
@click.option('--mesh_center_rot', required=True, type=float, help="the mesh rotation angle. (anti-clockwise)")
@click.option('--xis', required=True, type=str, help="xi0,xi1,dangle easting")
@click.option('--etas', required=True, type=str, help="eta0,eta1,dangle northing")
@click.option('--depths', required=True, type=str, help="depth0,depth1,ddepth")
@click.option('--peak_value', required=True, type=float, help="the peak value")
@click.option('--radius', required=True, type=float, help="the sphere radius in km")
@click.option('--output_path', required=True, type=str, help="the output list path")
def main(mesh_center_lat, mesh_center_lon, mesh_center_rot, xis, etas, depths, peak_value, radius, output_path):
    # * when zero rotation angle: xi -> easting, eta -> northing
    xi0, xi1, dangle = map(float, xis.split(","))
    eta0, eta1, dangle = map(float, etas.split(","))
    depth0, depth1, ddepth = map(float, depths.split(","))
    R_EARTH_KM = 6371.0

    # * local radial unit vector
    # note the mesh of specfem should set ellip to False
    x, y, z = lld2xyz(mesh_center_lon, mesh_center_lat, 0.0)
    vr0 = np.array([x, y, z])
    vr0 /= sum(vr0**2)**0.5

    # * local northing, easting unit vectors
    lat0 = np.arcsin(vr0[2])
    lon0 = np.arctan2(vr0[1], vr0[0])
    ve0 = np.array([-1*np.sin(lon0), np.cos(lon0), 0.0])
    vn0 = np.array([-1*np.sin(lat0)*np.cos(lon0), -1 *
                    np.sin(lat0)*np.sin(lon0), np.cos(lat0)])

    # * local basis along and perpendicular to the rotation azimuth
    az0 = -1*np.deg2rad(mesh_center_rot)
    va0 = vn0*np.cos(az0) + ve0*np.sin(az0)
    vp0 = -1*vn0*np.sin(az0) + ve0*np.cos(az0)

    # * ====== get grid of xi, eta: xi along va0, eta along vp0
    xi_list = np.arange(xi0, xi1+dangle/2, dangle)
    eta_list = np.arange(eta0, eta1+dangle/2, dangle)
    depth_list = np.arange(depth0, depth1+ddepth/2, ddepth)

    nxi = len(xi_list)
    neta = len(eta_list)
    ndepth = len(depth_list)

    npts = nxi*neta*ndepth
    xyz = np.zeros((3, npts))
    model_value = np.zeros(npts)

    dangle = np.deg2rad(dangle)

    i = 0
    for xi in xi_list:
        xi = np.deg2rad(xi)
        for eta in eta_list:
            eta = np.deg2rad(eta)
            v1 = np.cos(xi)*vp0 - np.sin(xi)*vr0
            v2 = np.cos(eta)*va0 - np.sin(eta)*vr0
            vx = np.cross(v1, v2)
            vx = vx/sum(vx**2)**0.5
            for depth in depth_list:
                sign = np.cos(np.pi*xi/dangle)*np.cos(np.pi*eta/dangle) * \
                    np.cos(np.pi*(depth-depth0)/ddepth)
                r = 1.0 - depth/R_EARTH_KM
                xyz[:, i] = r*vx
                model_value[i] = sign*peak_value
                i += 1

    # * ====== write out grid file
    with open(output_path, "w") as f:
        for i in range(npts):
            # f.write("%+14.7E  %+14.7E  %+14.7E  %.0f  %+8.5f\n" %
            #         (xyz[0, i], xyz[1, i], xyz[2, i], model_value[i]))
            f.write(
                f"{xyz[0,i]:+14.7E} {xyz[1,i]:+14.7E} {xyz[2,i]:+14.7E} {radius:.0f} {model_value[i]:+8.5f} \n")


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
