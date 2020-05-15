"""
Generate the string for making the gaussian mask.
"""
from glob import glob
from os.path import join, basename


def make_gaussian_mask(base_dir, sem_bin_dir, output_dir, nproc):
    all_event_dirs = sorted(glob(join(base_dir, "*")))
    command = join(sem_bin_dir, "xsem_make_gaussian_mask")
    result = ""
    result += "module load netcdf; "
    for index, each_event_dir in enumerate(all_event_dirs):
        gcmtid = basename(each_event_dir)
        mesh_dir = join(each_event_dir, "DATABASES_MPI")
        source_xyz_list = join(each_event_dir, "OUTPUT_FILES", "mask.xyz")
        out_dir = join(output_dir, gcmtid)
        result += f"mkdir -p {out_dir} ; date; "
        out_name = "mask_source"
        result += f"echo 'mask source for iter{index}:{gcmtid}' ; date; ibrun {command} {nproc} {mesh_dir} {source_xyz_list} {out_dir} {out_name}; "
    return result
