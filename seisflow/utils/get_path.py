from os.path import join


def get_asdf_fnames(gcmtid, min_periods, max_periods, data_asdf_directory, sync_asdf_directory):
    """
    get_asdf_fnames: get asdf fnames, min_periods=min_body,min_surface
    """
    # we have to make sure the asdf files in asdf_directory is complete.
    body_min_period, surface_min_period = min_periods.split(",")
    body_max_period, surface_max_period = max_periods.split(",")
    data_asdf_body_path = join(
        data_asdf_directory, f"{gcmtid}.preprocessed_{body_min_period}s_to_{body_max_period}s")
    data_asdf_surface_path = join(
        data_asdf_directory, f"{gcmtid}.preprocessed_{surface_min_period}s_to_{surface_max_period}s")
    sync_asdf_body_path = join(
        sync_asdf_directory, f"{gcmtid}.preprocessed_{body_min_period}s_to_{body_max_period}s")
    sync_asdf_surface_path = join(
        sync_asdf_directory, f"{gcmtid}.preprocessed_{surface_min_period}s_to_{surface_max_period}s")

    return data_asdf_body_path, sync_asdf_body_path, data_asdf_surface_path, sync_asdf_surface_path


def get_data_asdf_fnames(gcmtid, min_periods, max_periods, data_asdf_directory):
    """
    get_data_asdf_fnames: get asdf fnames, min_periods=min_body,min_surface
    """
    # we have to make sure the asdf files in asdf_directory is complete.
    body_min_period, surface_min_period = min_periods.split(",")
    body_max_period, surface_max_period = max_periods.split(",")
    data_asdf_body_path = join(
        data_asdf_directory, f"{gcmtid}.preprocessed_{body_min_period}s_to_{body_max_period}s")
    data_asdf_surface_path = join(
        data_asdf_directory, f"{gcmtid}.preprocessed_{surface_min_period}s_to_{surface_max_period}s")
    return data_asdf_body_path, data_asdf_surface_path
