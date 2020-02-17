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


def run_generate_model_perturbation_kernel(index, tags, database_list, base_directory):
    """
    generate the model perturbation. (the jobs should be run in parallel)
    """
    # ! notice the ending of the function is &
    julia_path = get_julia("scripts/get_perturbation.jl")
    result = ""
    # * init parameters
    nproc = int(database_list[index, 4])
    name = database_list[index, 0]
    target_basedir = join(base_directory, "absolute", name)
    reference_basedir = join(base_directory, "reference", name)
    mesh_basedir = join(base_directory, "mesh", name)
    output_basedir = join(base_directory, "perturbation", name)
    # * set up the command
    result += f"julia '{julia_path}' --target_basedir {target_basedir} --reference_basedir {reference_basedir} --mesh_basedir {mesh_basedir} \
        --output_basedir {output_basedir} --tags {tags} --nproc {nproc} & \n"
    return result


def run_interpolation(index, tags, database_list, base_directory):
    """
    run the interpolation command for the perturbed model, old: index, new: index+1.
    """
    # ! notice the ending of the function is ;
    julia_path = get_julia("specfem_gll.jl/src/program/xsem_interp_mesh2.jl")
    model_tags = tags
    result = ""
    # * init parameters
    database_list_old = database_list[index]
    database_list_new = database_list[index+1]
    nproc_old = int(database_list_old[4])
    nproc_new = int(database_list_new[4])
    old_mesh_dir = database_list_old[3]
    new_mesh_dir = database_list_new[3]
    # for the model dir, we should ref to the perturbation
    old_name = database_list_old[0]
    new_name = database_list_new[0]
    # notice, our old model dir would be the interpolated one if possible
    if(index == 0):
        old_model_dir = join(base_directory, "perturbation", old_name)
    else:
        past_name = database_list[index - 1, 0]
        old_model_dir = join(base_directory, "interpolation",
                             f"{past_name}__and__{old_name}")
    new_model_dir = join(base_directory, "perturbation", new_name)
    output_dir = join(base_directory, "interpolation",
                      f"{old_name}__and__{new_name}")
    # * now we construct the command
    result += f"ibrun julia '{julia_path}' --nproc_old {nproc_old} --old_mesh_dir {old_mesh_dir} --old_model_dir {old_model_dir} --nproc_new {nproc_new} --new_mesh_dir {new_mesh_dir} \
        --new_model_dir {new_model_dir} --model_tags {model_tags} --output_dir {output_dir}; \n"
    return result


def run_retrive_model(index, tags, database_list, base_directory):
    """
    Retrive the real final model from the model perturbation.
    """
    # target basedir is the last interpolated model
    target_basedir = join(base_directory, "interpolation",
                          f"{database_list[-2,0]}__and__{database_list[-1,0]}")
    # we have to keep all the reference_basedir the same (or cases like correction), use the first reference
    reference_basedir = database_list[0, 2]
