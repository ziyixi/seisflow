"""
Generate the string for making the gaussian mask.
"""
from glob import glob
from os.path import join


def make_gaussian_mask(base_dir, sem_bin_dir, nproc):
    all_event_dirs = glob(join(base_dir, "*"))
    command = join(sem_bin_dir, "xsem_make_gaussian_mask")
    result = ""
    result += "module load netcdf; "
    for each_event_dir in all_event_dirs:
        mesh_dir = join(each_event_dir, "DATABASES_MPI")
        source_xyz_list = join(each_event_dir, "OUTPUT_FILES", "mask.xyz")
        out_dir = join(each_event_dir, "DATABASES_MPI")
        out_name = "mask_source"
        result += f"ibrun {command} {nproc} {mesh_dir} {source_xyz_list} {out_dir} {out_name}; "
    return result
