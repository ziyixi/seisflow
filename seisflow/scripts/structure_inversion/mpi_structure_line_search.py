"""
mpi_structure_line_search.py: perform the structure line search in parallel for all the events.
"""
from mpi4py import MPI

from ...tasks.structure_inversion import line_search
from ...utils.asdf_io import VirAsdf
from ...utils.setting import SURFACE_THRESHOLD
from ...utils.get_path import get_data_asdf_fnames, get_sync_asdf_fnames
from ...utils.load_files import load_windows
from glob import glob
from os.path import join, basename
import numpy as np

comm = MPI.COMM_WORLD  # pylint: disable=c-extension-no-member
size = comm.Get_size()
rank = comm.Get_rank()


def get_gcmtid_this_rank(cmts_directory):
    """
    Get gcmt ids used in this rank.
    """
    all_gcmtid_files = sorted(glob(join(cmts_directory, "*")))
    all_gcmtids = [basename(item) for item in all_gcmtid_files]
    assert len(all_gcmtids) == size
    return all_gcmtids[rank]


def prepare_line_search_parameters(gcmtid_this_rank, data_asdf_directory, sync_asdf_raw_directory, sync_asdf_perturbed_directory,
                                   windows_directory, data_info_directory, stations_path, min_periods, max_periods):
    """
    Prepare parameters for the structure line search.
    """
    # * firstly we get paths for data_virasdf_body, data_virasdf_surface, sync_virasdf_body_raw,
    # * sync_virasdf_surface_raw, sync_virasdf_body_perturbed, sync_virasdf_surface_perturbed
    data_virasdf_body_path, data_virasdf_surface_path = get_data_asdf_fnames(
        gcmtid_this_rank, min_periods, max_periods, data_asdf_directory)
    sync_virasdf_body_raw_path, sync_virasdf_surface_raw_path = get_sync_asdf_fnames(
        gcmtid_this_rank, min_periods, max_periods, sync_asdf_raw_directory)
    sync_virasdf_body_perturbed_path, sync_virasdf_surface_perturbed_path = get_sync_asdf_fnames(
        gcmtid_this_rank, min_periods, max_periods, sync_asdf_perturbed_directory)
    # and now we can load asdf files.
    data_virasdf_body = VirAsdf()
    data_virasdf_body.read_asdf(data_virasdf_body_path)
    data_virasdf_surface = VirAsdf()
    data_virasdf_surface.read_asdf(data_virasdf_surface_path)
    sync_virasdf_body_raw = VirAsdf()
    sync_virasdf_body_raw.read_asdf(sync_virasdf_body_raw_path)
    sync_virasdf_surface_raw = VirAsdf()
    sync_virasdf_surface_raw.read_asdf(sync_virasdf_surface_raw_path)
    sync_virasdf_body_perturbed = VirAsdf()
    sync_virasdf_body_perturbed.read_asdf(sync_virasdf_body_perturbed_path)
    sync_virasdf_surface_perturbed = VirAsdf()
    sync_virasdf_surface_perturbed.read_asdf(
        sync_virasdf_surface_perturbed_path)
    # * we can load the windows
    windows = load_windows(gcmtid_this_rank, windows_directory)
    # * we should get if consider_surface
    # firstly we should get the value of the depth
    depth = data_virasdf_body.get_events()[0].preferred_origin().depth / 1000.0
    if (depth <= SURFACE_THRESHOLD):
        consider_surface = True
    else:
        consider_surface = False
