"""
interp_new_model.py: combine several models as a new hybrid model.
"""
import numpy as np
from ... import get_julia
import sh
from os.path import join


def load_databases_list(database_list_path):
    """
    load databases list: name absolute_dir ref_dir mesh_dir nproc (in each row)
    """
    database_list = np.loadtxt(database_list_path, dtype=np.str)
    return database_list


def secure_link(from_path, to_path):
    """
    secure link to prevent the case if the to_path exists.
    """
    try:
        sh.unlink(to_path)
    except sh.ErrorReturnCode_1:
        pass
    sh.ln("-s", from_path, to_path)


def init_structure(base_directory, database_list):
    """
    init the structure for generating the model.
    """
    # * make some base directories
    sh.mkdir("-p", base_directory)
    sh.mkdir("-p", join(base_directory, "absolute"))
    sh.mkdir("-p", join(base_directory, "reference"))
    sh.mkdir("-p", join(base_directory, "mesh"))
    sh.mkdir("-p", join(base_directory, "perturbation"))
    sh.mkdir("-p", join(base_directory, "interpolation"))
    # the generated directories will contain the final model
    sh.mkdir("-p", join(base_directory, "generated"))
    # * for each row in database_list, make soft links
    for each_row in database_list:
        name = each_row[0]
        absolute_path = each_row[1]
        reference_path = each_row[2]
        mesh_path = each_row[3]
        # make soft links
        secure_link(absolute_path, join(base_directory, "absolute", name))
        secure_link(reference_path, join(base_directory, "reference", name))
        secure_link(mesh_path, join(base_directory, "mesh", name))
        # make perturbation directories
        sh.mkdir("-p", join(base_directory, "perturbation", name))
    # * generate the interpolation directories
    n_model = database_list.shape[0]
    for imodel in range(n_model - 1):
        name_this = database_list[imodel, 0]
        name_next = database_list[imodel, 1]
        sh.mkdir("-p", join(base_directory, "interpolation",
                            f"{name_this}__and__{name_next}"))


def run_generate_model_perturbation(name, nproc, tags, base_directory):
    """
    generate the model perturbation. (the jobs should be run in parallel)
    """
    julia_path = get_julia("scripts/get_perturbation.jl")
    result = ""
    # * init parameters
    nproc = int(nproc)
    target_basedir = join(base_directory, "absolute", name)
    reference_basedir = join(base_directory, "reference", name)
    mesh_basedir = join(base_directory, "mesh", name)
    output_basedir = join(base_directory, "perturbation", name)
    result += f"julia {julia_path} --target_basedir {target_basedir} --reference_basedir {reference_basedir} --mesh_basedir {mesh_basedir} \
        --output_basedir {output_basedir} --tags {tags} --nproc {nproc} & \n"
    return result


def run_interpolation(index, database_list, base_directory):
    """
    run the interpolation command for the perturbed model.
    """
    julia_path = get_julia("scripts/get_perturbation.jl")
