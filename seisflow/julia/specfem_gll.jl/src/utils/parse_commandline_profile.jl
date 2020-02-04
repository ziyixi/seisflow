using ArgParse

function parse_commandline()
    s = ArgParseSettings(description = """This program is used to generate profile lines""")
    @add_arg_table s begin
        """--nproc_old"""
        help = "number of slices of the old mesh"
        arg_type = Int64
        required = true
        """--old_mesh_dir"""
        help = "directory holds proc*_reg1_solver_data.bin"
        arg_type = String
        required = true
        """--old_model_dir"""
        help = "directory holds proc*_reg1_<model_tag>.bin"
        arg_type = String
        required = true
        """--model_tags"""
        help = "comma delimited string, e.g. vsv,vsh,rho"
        arg_type = String
        required = true
        """--output_file"""
        help = "output directory for interpolated model files"
        arg_type = String
        required = true
        """--region"""
        help = "lon1/lat1/lon2/lat2/dep1/dep2, if dep1==dep2, cut horizontally, otherwise cut vertically"
        arg_type = String
        required = true
        """--npts"""
        help = "hnpts/vnpts if cut vertically, lon_npts/lat_npts if cut horizontally, should be interger"
        arg_type = String
        required = true
    end
    return parse_args(s)
end
