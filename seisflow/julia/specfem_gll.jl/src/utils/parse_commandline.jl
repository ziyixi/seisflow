using ArgParse

function parse_commandline()
    s = ArgParseSettings(description = """1) This program must run in parallel, e.g. mpirun -n <nproc> ...
    2) use values_from_mesher.h for the old mesh to compile
    use two values: ANGULAR_WIDTH_XI_IN_DEGREES_VAL, NEX_XI_VAL
    3) the new_model is used as the background model where the
    old mesh doesn't cover""")
    @add_arg_table s begin
        "--nproc_old"
        help = "number of slices of the old mesh"
        arg_type = Int64
        required = true
        "--old_mesh_dir"
        help = "directory holds proc*_reg1_solver_data.bin"
        arg_type = String
        required = true
        "--old_model_dir"
        help = "directory holds proc*_reg1_<model_tag>.bin"
        arg_type = String
        required = true
        "--nproc_new"
        help = "number of slices of the new mesh"
        arg_type = Int64
        required = true
        "--new_mesh_dir"
        help = "directory holds proc*_reg1_solver_data.bin"
        arg_type = String
        required = true
        "--new_model_dir"
        help = "directory for new model files as background model"
        arg_type = String
        required = true
        "--model_tags"
        help = "comma delimited string, e.g. vsv,vsh,rho"
        arg_type = String
        required = true
        "--output_dir"
        help = "output directory for interpolated model files"
        arg_type = String
        required = true
    end
    return parse_args(s)
end
