"""
Convert the processed asdf file to the sac format.
"""
from glob import glob
from os.path import basename, join
import pyasdf
import obspy
import pickle
import numpy as np
import click


def load_pickle(pickle_path):
    with open(pickle_path, "rb") as f:
        data = pickle.load(f)
    return data


@click.command()
@click.option('--asdf_path', required=True, type=str)
@click.option('--data_info_dir', required=True, type=str)
@click.option('--stations_fname', required=True, type=str)
@click.option('--output_dir', required=True, type=str)
@click.option('--az_range', required=True, type=str)
def main(asdf_path, data_info_dir, stations_fname, output_dir, az_range):
    # load stations info
    stations = np.loadtxt(stations_fname, dtype=np.str)
    stations_mapper = {}
    for row in stations:
        net_sta = f"{row[1]}.{row[0]}"
        stations_mapper[net_sta] = (
            float(row[2]), float(row[3]), float(row[4]))
    # load header info
    az_dict = load_pickle(join(data_info_dir, "extra.az.pkl"))
    baz_dict = load_pickle(join(data_info_dir, "extra.baz.pkl"))
    evdp_dict = load_pickle(join(data_info_dir, "extra.evdp.pkl"))
    evla_dict = load_pickle(join(data_info_dir, "extra.evla.pkl"))
    evlo_dict = load_pickle(join(data_info_dir, "extra.evlo.pkl"))
    gcarc_dict = load_pickle(join(data_info_dir, "extra.gcarc.pkl"))
    t1_dict = load_pickle(join(data_info_dir, "traveltime.P.pkl"))
    t2_dict = load_pickle(join(data_info_dir, "traveltime.S.pkl"))

    asdf_file = pyasdf.ASDFDataSet(asdf_path, mode="r")
    thebasename = basename(asdf_path)
    gcmtid = thebasename.split(".")[0]
    tag = thebasename.split(".")[1]
    avaliable_net_sta = set(list(stations_mapper.keys())) & set(
        asdf_file.waveforms.list())
    az_range = map(float, az_range.split(","))
    for net_sta in sorted(avaliable_net_sta):
        if(az_range[0] <= az_dict[gcmtid][net_sta] <= az_range[1]):
            # z
            tr = asdf_file.waveforms[net_sta][tag].select(component="Z")[
                0].copy()
            tr.write(join(output_dir, tr.id), format="SAC")
            tr = obspy.read(join(output_dir, tr.id))[0]
            tr.stats.sac.az = az_dict[gcmtid][net_sta]
            tr.stats.sac.baz = baz_dict[gcmtid][net_sta]
            tr.stats.sac.evdp = evdp_dict[gcmtid][net_sta]
            tr.stats.sac.evla = evla_dict[gcmtid][net_sta]
            tr.stats.sac.evlo = evlo_dict[gcmtid][net_sta]
            tr.stats.sac.gcarc = gcarc_dict[gcmtid][net_sta]
            tr.stats.sac.t1 = t1_dict[gcmtid][net_sta]
            tr.stats.sac.t2 = t2_dict[gcmtid][net_sta]
            tr.stats.sac.kt1 = "P"
            tr.stats.sac.kt2 = "S"
            tr.stats.sac.stla = stations_mapper[net_sta][0]
            tr.stats.sac.stlo = stations_mapper[net_sta][1]
            tr.stats.sac.stel = stations_mapper[net_sta][2]
            tr.stats.sac.stdp = 0.0
            tr.write(join(output_dir, tr.id), format="SAC")

            # r
            tr = asdf_file.waveforms[net_sta][tag].select(component="R")[
                0].copy()
            tr.write(join(output_dir, tr.id), format="SAC")
            tr = obspy.read(join(output_dir, tr.id))[0]
            tr.stats.sac.az = az_dict[gcmtid][net_sta]
            tr.stats.sac.baz = baz_dict[gcmtid][net_sta]
            tr.stats.sac.evdp = evdp_dict[gcmtid][net_sta]
            tr.stats.sac.evla = evla_dict[gcmtid][net_sta]
            tr.stats.sac.evlo = evlo_dict[gcmtid][net_sta]
            tr.stats.sac.gcarc = gcarc_dict[gcmtid][net_sta]
            tr.stats.sac.t1 = t1_dict[gcmtid][net_sta]
            tr.stats.sac.t2 = t2_dict[gcmtid][net_sta]
            tr.stats.sac.kt1 = "P"
            tr.stats.sac.kt2 = "S"
            tr.stats.sac.stla = stations_mapper[net_sta][0]
            tr.stats.sac.stlo = stations_mapper[net_sta][1]
            tr.stats.sac.stel = stations_mapper[net_sta][2]
            tr.stats.sac.stdp = 0.0
            tr.write(join(output_dir, tr.id), format="SAC")

            # t
            tr = asdf_file.waveforms[net_sta][tag].select(component="T")[
                0].copy()
            tr.write(join(output_dir, tr.id), format="SAC")
            tr = obspy.read(join(output_dir, tr.id))[0]
            tr.stats.sac.az = az_dict[gcmtid][net_sta]
            tr.stats.sac.baz = baz_dict[gcmtid][net_sta]
            tr.stats.sac.evdp = evdp_dict[gcmtid][net_sta]
            tr.stats.sac.evla = evla_dict[gcmtid][net_sta]
            tr.stats.sac.evlo = evlo_dict[gcmtid][net_sta]
            tr.stats.sac.gcarc = gcarc_dict[gcmtid][net_sta]
            tr.stats.sac.t1 = t1_dict[gcmtid][net_sta]
            tr.stats.sac.t2 = t2_dict[gcmtid][net_sta]
            tr.stats.sac.kt1 = "P"
            tr.stats.sac.kt2 = "S"
            tr.stats.sac.stla = stations_mapper[net_sta][0]
            tr.stats.sac.stlo = stations_mapper[net_sta][1]
            tr.stats.sac.stel = stations_mapper[net_sta][2]
            tr.stats.sac.stdp = 0.0
            tr.write(join(output_dir, tr.id), format="SAC")

            # # PZ file
            # stationxml = asdf_file.waveforms[net_sta].StationXML
            # stationxml.write(join(output_dir, f"{net_sta}.PZ"), format="SACPZ")


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
