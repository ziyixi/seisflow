import numpy as np
import click


def convert(x, y, z):
    R = 6371
    dep = R - R * np.sqrt(x ** 2 + y ** 2 + z ** 2)
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    rh = np.sqrt(x ** 2 + y ** 2)
    theta = np.rad2deg(np.arcsin(rh / r))
    lat = 90 - theta
    phi = np.rad2deg(np.arctan(x / y))
    lon = 90-phi
    return lon, lat, dep


@click.command()
@click.option('--psf_list_path', required=True, type=str)
@click.option('--output_path', required=True, type=str)
def main(psf_list_path, output_path):
    psf_list = np.loadtxt(psf_list_path)
    with open(output_path, "w") as f:
        for row in psf_list:
            lon, lat, dep = convert(row[0], row[1], row[2])
            f.write(f"{lon:.2f} {lat:.2f} {dep:.0f} \n")


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
