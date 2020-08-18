using ArgParse
using ProgressMeter

include("../utils/readfiles.jl")

function parse_commandline()
    s = ArgParseSettings(description="read through d670 kernels and write to a text file with each row as: lon, lat, dep, kernel_val")
    @add_arg_table! s begin
        "--nproc"
        help = "the number of slices in the mesh"
        arg_type = Int64
        required = true
        "--base_dir"
        help = "the DATABASE directory"
        arg_type = String
        required = true
        "--output_path"
        help = "the output text file path"
        arg_type = String
        required = true     
    end   
    return parse_args(s)
end


function xyz2lld(x::Float64, y::Float64, z::Float64)
    R = 6371.0
    r = sqrt(x^2 + y^2 + z^2)
    dep = R - R * r
    rh = sqrt(x^2 + y^2)
    θ = asind(rh / r)
    lat = 90 - θ
    ϕ = atand(x / y)
    lon = 90 - ϕ
    return [lon,lat,dep]
end


function write2file(output_path, kernel_d670_all, position_d670_all, NSPEC2D_670::Int64, nproc::Int64)
    open(output_path, "w") do io
        for iproc in 1:nproc
            for ispec in 1:NSPEC2D_670
                for iglly in 1:NGLLY
                    for igllx in 1:NGLLX
                        kernel_val = kernel_d670_all[igllx,iglly,ispec,iproc]
                        lon, lat, dep = position_d670_all[:,igllx,iglly,ispec,iproc]
                        write(io, "$lon $lat $dep $kernel_val \n")
                    end
                end
            end
        end
end
end

function main()
    # parse command line
    command_args = parse_commandline()  
    nproc = command_args["nproc"]
    base_dir = command_args["base_dir"]
    output_path = command_args["output_path"]

    # read d670 kernel and boundary mesh
    # for iproc, read the mesh file to locate the location
    boudary_disc_data = sem_boundary_disc_read(base_dir, 0)
    kernel_d670_all = zeros(Float64, NGLLX, NGLLY, boudary_disc_data.NSPEC2D_670, nproc)
    position_d670_all = zeros(Float64, 3, NGLLX, NGLLY, boudary_disc_data.NSPEC2D_670, nproc)
    @showprogress for iproc in 1:nproc
        kernel_d670_all[:,:,:,iproc] = sem_d670_read(base_dir, iproc - 1, boudary_disc_data.NSPEC2D_670)
        mesh_data_iproc = sem_mesh_data()
        mesh_data_iproc = sem_mesh_read(base_dir, iproc - 1)
        for ispec in ibelm_670_top
            for igllx in 1:NGLLX
                for iglly in 1:NGLLY
                    iglob = mesh_old.ibool[igllx,iglly,1,ispec]
                    x, y, z = mesh_data_iproc.xyz_glob[:,iglob]
                    position_d670_all[:,igllx,iglly,ispec,iproc] = xyz2lld(x, y, z)
                end
            end
        end
    end

    # now we write to the file
    write2file(output_path, kernel_d670_all, position_d670_all, boudary_disc_data.NSPEC2D_670, nproc)
end

main()