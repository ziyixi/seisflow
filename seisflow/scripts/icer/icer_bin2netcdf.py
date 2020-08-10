"""
bin2netcdf.py: submit a job to convert the mesh bin files to the netcdf file.
"""

import sys
from os.path import basename, dirname, join, expanduser

import click
import sh

from ... import get_julia
from ...slurm.submit_job import submit_job


def bin2ppm(nproc_old, model_tags, region, npts, nproc,
            old_mesh_dir, old_model_dir, output_dir):
    """
    convert the bin files to the ppm model.
    """
    result = ""
    julia_path = get_julia("specfem_gll.jl/src/program/get_ppm_model.jl")
    latnproc, lonnproc = map(int, nproc.split("/"))
    nproc_ppm2netcdf = latnproc * lonnproc
    # ! note there is a issue of precompiling the code in a race condition, refer to https://github.com/simonbyrne/PkgLock.jl to solve the problem
    # result += "julia --project -e 'push!(LOAD_PATH, \"@pkglock\"); using PkgLock; PkgLock.instantiate_precompile()'\n"
    result += "module purge;module load GCC/8.2.0-2.31.1;module load OpenMPI/3.1.3;"
    result += f"srun -n {nproc_ppm2netcdf} julia '{julia_path}' --nproc_old {nproc_old} --old_mesh_dir {old_mesh_dir} --old_model_dir {old_model_dir} --model_tags {model_tags} --output_file {output_dir} --region {region} --npts {npts} --nproc {nproc}; \n"
    return result


def ppm2netcdf(base_dir, region, npts, nproc, parameters, out_path, history):
    """
    convert the ppm model to the netcdf model.
    """
    result = ""
    pyexec = sys.executable
    result += f"srun -n {nproc} {pyexec} -m seisflow.scripts.mesh.convert_txt_output_2_netcdf --base_dir {base_dir} --region {region} --npts {npts} --nproc {nproc} \
        --parameters {parameters} --out_path {out_path} --history '{history}'; \n"
    return result


@click.command()
@click.option('--nproc_old', required=True, type=int, help="number of processes used in the bin files")
@click.option('--model_tags', required=True, type=str, help="the model tags, eg: vpv,vph")
@click.option('--region', required=True, type=str, help="the region to interpolate, lon1/lat1/lon2/lat2/dep1/dep2")
@click.option('--npts', required=True, type=str, help="nlon/nlat/ndep")
@click.option('--nproc', required=True, type=str, help="latnproc/lonnproc")
@click.option('--history', required=True, type=str, help="the info to write into the netcdf file")
@click.option('--time', required=True, type=str, help="the time to run as a job")
@click.option('--old_mesh_dir', required=True, type=str, help="the mesh directory")
@click.option('--old_model_dir', required=True, type=str, help="the bin file model directory")
@click.option('--output_path', required=True, type=str, help="the output netcdf path")
def main(nproc_old, model_tags, region, npts, nproc, history, time,
         old_mesh_dir, old_model_dir, output_path):
    """
    submit a job to convert the mesh bin files to the netcdf file.
    """
    # * firstly we have to make a tmp file to store the
    # temp_directory = tempfile.mkdtemp()
    temp_directory = join(dirname(output_path), "."+basename(output_path))
    sh.mkdir("-p", temp_directory)
    # * generate the ppm model
    result = ""
    # ! fix MPI cache issue
    home = expanduser("~")
    julia_path = join(home, ".julia")
    result += f"export JULIA_DEPOT_PATH={julia_path}\n"
    result += f"TMPDIR=`mktemp -d`\n"
    result += f"mkdir $TMPDIR/compiled\n"
    # if use v1.1
    result += f"rsync -au $JULIA_DEPOT_PATH/compiled/v1.1 $TMPDIR/compiled/\n"
    result += f"export JULIA_DEPOT_PATH=$TMPDIR:$JULIA_DEPOT_PATH\n"
    result += "date; \n"
    # bin2ppm npts use nlat/nlon/ndep
    nlon, nlat, ndep = npts.split("/")
    npts_bin2ppm = "/".join([nlat, nlon, ndep])
    result += bin2ppm(nproc_old, model_tags, region, npts_bin2ppm, nproc,
                      old_mesh_dir, old_model_dir, temp_directory)
    # * convert to netcdf model
    # region in bin2ppm: lon1/lat1/lon2/lat2/dep1/dep2; region in ppm2netcdf: minlon/maxlon/minlat/maxlat/mindep/maxdep
    # we can treat 1 as the min, 2 as the max
    result += "date; \n"
    minlon, minlat, maxlon, maxlat, mindep, maxdep = region.split("/")
    region_ppm2netcdf = "/".join([minlon, maxlon,
                                  minlat, maxlat, mindep, maxdep])
    # here the value of nproc should be latnproc*lonnproc
    latnproc, lonnproc = map(int, nproc.split("/"))
    nproc_ppm2netcdf = latnproc * lonnproc
    # ! we should change npts from nlon/nlat/ndep to latnpts/lonnpts/vnpts
    result += ppm2netcdf(temp_directory, region_ppm2netcdf, npts, nproc_ppm2netcdf,
                         model_tags, output_path, history)
    result += "date; \n"
    # * submit the job
    submit_job("bin2netcdf", result, None, nproc_ppm2netcdf,
               None, time, None, "icer_flexiable")


if __name__ == "__main__":
    main()  # pylint: disable=no-value-for-parameter
