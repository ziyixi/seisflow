"""
interp_new_model.py: combine several models as a new hybrid model.
"""
from os.path import join, expanduser

import click
import numpy as np
import sh

from ... import get_julia
from ...slurm.submit_job import submit_job


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
    for imodel in range(1, n_model):
        names = database_list[:(imodel+1), 0]
        sh.mkdir("-p", join(base_directory, "interpolation",
                            "__and__".join(names)))


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
    result += f"ibrun -n 1 julia -t 68 '{julia_path}' --target_basedir {target_basedir} --reference_basedir {reference_basedir} --mesh_basedir {mesh_basedir} \
        --output_basedir {output_basedir} --tags {tags} --nproc {nproc} ; \n"
    return result


def run_interpolation(index, tags, database_list, base_directory):
    """
    run the interpolation command for the perturbed model, old: index+1, new: index.
    """
    # ! notice the ending of the function is ;
    julia_path = get_julia("specfem_gll.jl/src/program/xsem_interp_mesh2.jl")
    result = ""
    # * we can get the values used in the command
    # * new
    nproc_new = int(database_list[0, 4])
    new_mesh_dir = join(base_directory, "mesh", database_list[0, 0])
    if (index == 0):
        new_model_dir = join(
            base_directory, "perturbation", database_list[index, 0])
    else:
        names = database_list[: (index + 1), 0]
        new_model_dir = join(base_directory, "interpolation",
                             "__and__".join(names))
    # * old
    nproc_old = int(database_list[index + 1, 4])
    old_mesh_dir = join(base_directory, "mesh", database_list[index + 1, 0])
    old_model_dir = join(base_directory, "perturbation",
                         database_list[index + 1, 0])
    # * now we construct the command
    output_names = database_list[: (index + 2), 0]
    output_dir = join(base_directory, "interpolation",
                      "__and__".join(output_names))
    result += f"ibrun julia '{julia_path}' --nproc_old {nproc_old} --old_mesh_dir {old_mesh_dir} --old_model_dir {old_model_dir} --nproc_new {nproc_new} --new_mesh_dir {new_mesh_dir} \
        --new_model_dir {new_model_dir} --model_tags {tags} --output_dir {output_dir}; \n"
    return result


def run_retrive_model(tags, database_list, base_directory):
    """
    Retrive the real final model from the model perturbation.
    """
    # target basedir is the last interpolated model
    target_basedir = join(base_directory, "interpolation",
                          "__and__".join(database_list[:, 0]))
    # we have to keep all the reference_basedir the same (or cases like correction), use the first reference
    reference_basedir = join(base_directory, "reference", database_list[0, 0])
    mesh_basedir = join(base_directory, "mesh", database_list[0, 0])
    output_basedir = join(base_directory, "generated")
    nproc = int(database_list[0, 4])
    # * build up the command
    julia_path = get_julia("scripts/retrive_model.jl")
    result = ""
    result += f"ibrun -n 1 julia -t 68 '{julia_path}' --target_basedir {target_basedir} --reference_basedir {reference_basedir} --mesh_basedir {mesh_basedir} --output_basedir {output_basedir}  \
        --tags {tags} --nproc {nproc};\n "
    return result


@click.command()
@click.option('--base_directory', required=True, type=str, help="the base interpolation directory")
@click.option('--database_list_path', required=True, type=str, help="the database_list file path")
@click.option('--tags', required=True, type=str, help="the tags to interpolate, eg: vsv,vsh")
@click.option('--n_node', required=True, type=int, help="the number of nodes to use")
@click.option('--ntasks', required=True, type=int, help="the number of tasks to use")
@click.option('--partition', required=True, type=str, help="the partition to use, eg: skx-normal")
@click.option('--time', required=True, type=str, help="the job time")
@click.option('--account', required=True, type=str, help="the stampede2 account")
def main(base_directory, database_list_path, tags, n_node, ntasks,
         partition, time, account):
    """
    Build up the hybrid model depends on previous gll models.
    """
    database_list = load_databases_list(database_list_path)
    # ! note we should reset the values in constants.jl
    result = ""
    # ! fix MPI cache issue
    home = expanduser("~")
    julia_path = join(home, ".julia")
    result += f"export JULIA_DEPOT_PATH={julia_path}\n"
    result += f"TMPDIR=`mktemp -d`\n"
    result += f"mkdir $TMPDIR/compiled\n"
    # if use v1.1
    result += f"rsync -au $JULIA_DEPOT_PATH/compiled/v1.6 $TMPDIR/compiled/\n"
    result += f"export JULIA_DEPOT_PATH=$TMPDIR:$JULIA_DEPOT_PATH\n"
    result += "date; \n"
    # * firstly we generate the structure of this workflow.
    init_structure(base_directory, database_list)
    # * then we generate the model perturbations.
    nmodel, _ = database_list.shape
    for imodel in range(nmodel):
        result += run_generate_model_perturbation_kernel(
            imodel, tags, database_list, base_directory)
    result += "date; \n"
    # * we interpolate the model in order.
    for imodel in range(nmodel - 1):
        # the nproc should be corresponding to the nproc of the first row.
        result += run_interpolation(imodel, tags,
                                    database_list, base_directory)
        result += "date; \n"
    # * now we retrive the model to the absolute values
    result += run_retrive_model(tags, database_list, base_directory)
    result += "date; \n"
    # * now we submit the job
    submit_job("interp", result, n_node, ntasks,
               partition, time, account, "stampede2")


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
