"""
get_paraview_slice_axis.py: calculate the axis used in slicing paraview moddels.
"""
import numpy as np
from sympy import Plane, Point3D


def sind(x):
    """
    sind: sin(x) in degree.
    """
    return np.sin(np.deg2rad(x))


def cosd(x):
    """
    cosd: cos(x) in degree.
    """
    return np.cos(np.deg2rad(x))


def convert_ll2xyz(latitude, longitude):
    """
    convert_ll2xyz: convert the latitude and longitude (degree) to x,t,z (normalized)
    """
    theta = 90 - latitude
    phi = longitude
    x = sind(theta) * cosd(phi)
    y = sind(theta) * sind(phi)
    z = cosd(theta)
    return x, y, z


if __name__ == "__main__":
    import click

    @click.command()
    @click.option('--p1', required=True, type=str, help="point1 lon1,lat1")
    @click.option('--p2', required=True, type=str, help="point2 lon2,lat2")
    def main(p1, p2):
        """
        The main function.
        """
        lon1, lat1 = map(float, p1.split(","))
        lon2, lat2 = map(float, p2.split(","))
        x1, y1, z1 = convert_ll2xyz(lat1, lon1)
        x2, y2, z2 = convert_ll2xyz(lat2, lon2)
        slice_plane = Plane(Point3D(x1, y1, z1), Point3D(
            x2, y2, z2), Point3D(0, 0, 0))
        axis = slice_plane.normal_vector
        # print out the result
        print(f"x: {float(axis[0])}")
        print(f"y: {float(axis[1])}")
        print(f"z: {float(axis[2])}")
        print(f"point1 {x1} {y1} {z1}")
        print(f"point2 {x2} {y2} {z2}")

    main()  # pylint: disable=no-value-for-parameter
