using ArgParse
include("../specfem_gll.jl/tasks/numerical_operation_on_bin_file.jl")
include("../specfem_gll.jl/src/utils/readfiles.jl")

function parse_command_line()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--target_basedir"
            help = "the target gll directory"
        "--reference_basedir"
            help = "the reference gll directory"
        "--mesh_basedir"
            help = "the mesh directory"
        "--output_basedir"
            help = "the output gll directory for perturbation"
        "--tags"
            help = "the tags, eg: alpha,beta"
        "--nproc"
            help = "number of processors the gll directory correspons to"
    end
    return parse_args(s)
    end

function main()
    # parse args
    parsed_args = parse_command_line()
    target_basedir = parsed_args["target_basedir"]
    mesh_basedir = parsed_args["mesh_basedir"]
    reference_basedir = parsed_args["reference_basedir"]
    output_basedir = parsed_args["output_basedir"]
    tags = parsed_args["tags"]
    nproc = parse(Int64, parsed_args["nproc"])
    # get nspec
    mesh_info = sem_mesh_read(mesh_basedir, 0)
    nspec = mesh_info.nspec
    # run generate_perturbation
    generate_perturbation(target_basedir, reference_basedir, output_basedir, tags, nproc, nspec)
end

main()