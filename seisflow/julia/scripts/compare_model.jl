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
        "--tag"
            help = "the tag"
        "--iproc"
            help = "the bin file index"
    end
    return parse_args(s)
    end

function main()
    # parse args
    parsed_args = parse_command_line()
    target_basedir = parsed_args["target_basedir"]
    mesh_basedir = parsed_args["mesh_basedir"]
    reference_basedir = parsed_args["reference_basedir"]
    tag = parsed_args["tag"]
    iproc = parse(Int64, parsed_args["iproc"])
    # get nspec
    mesh_info = sem_mesh_read(mesh_basedir, 0)
    nspec = mesh_info.nspec
    # run generate_perturbation
    compare_model(target_basedir::String, reference_basedir::String, tag::String,  nspec::Int64, iproc::Int64)
end

main()