import pickle
from os.path import join

import pyasdf


def save_pickle_event(to_save, output_dir, used_gcmtid):
    output_path = join(output_dir, f"{used_gcmtid}.pkl")
    with open(output_path, "wb") as f:
        pickle.dump(to_save, f)


def save_cmtsolution(output_path, cmtsolution):
    """
    save cmtsolution file.
    """
    cmtsolution.write(output_path, format="CMTSOLUTION")
    # remove the final line (or specfem will report error)
    with open(output_path, "r") as f:
        lines = f.readlines()
    with open(output_path, "w") as f:
        for i in range(12):
            f.write(lines[i])
        f.write(lines[-1].split("\n")[0])


def save_adjoint_to_asdf(adjoint_source_zerolagcc, output_directory, gcmtid, raw_sync_virasdf):
    """
    save_adjoint_to_asdf: save the adjoint source to asdf format to be read by Specfem.
    """
    # * note that the components name may not be like MX*
    waveforms = raw_sync_virasdf.get_waveforms()
    output_fname = join(output_directory, f"{gcmtid}.h5")
    components = ["E", "N", "Z"]
    with pyasdf.ASDFDataSet(output_fname, mode="w", mpi=False) as output_asdf:
        for net_sta in adjoint_source_zerolagcc:
            raw_st = waveforms[net_sta]["st"]
            head_component = raw_st[0].stats.channel[:2]
            for index_component in range(3):
                specfem_adj_source = adjoint_source_zerolagcc[net_sta][index_component, :]
                tag = net_sta.replace(".", "_") + "_" + \
                    head_component+components[index_component]
                output_asdf.add_auxiliary_data(
                    data=specfem_adj_source, data_type="AdjointSources", path=tag, parameters={})
