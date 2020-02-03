using ArgParse
include("../specfem_gll.jl/scripts/perturbation_bin_file.jl")
include("../specfem_gll.jl/src/utils/readfiles.jl")

function parse_command_line()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--target_basedir"
            help = "the target gll directory"
        "--reference_basedir"
            help = "the reference gll directory"
        "--output_basedir"
            help = "the output gll directory for the real model"
        "--nproc"
            help = "number of processors the gll directory correspons to"
    end
    return parse_args(s)
end

function main()
    # parse args
    parsed_args = parse_command_line()
    target_basedir = parsed_args["target_basedir"]
    reference_basedir = parsed_args["reference_basedir"]
    output_basedir = parsed_args["output_basedir"]
    nproc = parse(Int64, parsed_args["nproc"])
    # get nspec
    mesh_info = sem_mesh_read(reference_basedir, 0)
    nspec = mesh_info.nspec
    # run generate_perturbation
    generate_real(target_basedir, reference_basedir, output_basedir, nproc, nspec)
end

main()