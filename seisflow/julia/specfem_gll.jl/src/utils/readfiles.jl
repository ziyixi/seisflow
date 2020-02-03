using FortranFiles
include("../setting/constants.jl")
include("../types/types.jl")

# * mesh
# nspec for my simulation setting is 4480
function sem_mesh_read(basedir::String, iproc::Int64)
    # * init
    mesh_data = sem_mesh_data()
    f = sem_io_open_file_for_read(basedir, iproc, "solver_data") 

    nspec = Int64(read(f, Int32))
    nglob = Int64(read(f, Int32))

    sem_mesh_init!(mesh_data, nspec, nglob)
    mesh_data.nspec = Int64(nspec)
    mesh_data.nglob = Int64(nglob)

    dummy = zeros(Float32, nglob)
    dummy4 = zeros(Float32, NGLLX, NGLLY, NGLLZ, nspec)

    read(f, dummy)
    mesh_data.xyz_glob[1,:] .= Float64.(dummy)
    read(f, dummy)
    mesh_data.xyz_glob[2,:] .= Float64.(dummy)
    read(f, dummy)
    mesh_data.xyz_glob[3,:] .= Float64.(dummy)

    dummy_int1 = similar(mesh_data.ibool, Int32)
    read(f, dummy_int1)
    mesh_data.ibool = Int64.(dummy_int1)

    dummy_int2 = similar(mesh_data.idoubling, Int32)
    read(f, dummy_int2)
    mesh_data.idoubling = Int64.(dummy_int2)

    read(f, mesh_data.ispec_is_tiso)

    read(f, dummy4)
    mesh_data.xix .= Float64.(dummy4)
    read(f, dummy4)
    mesh_data.xiy .= Float64.(dummy4)
    read(f, dummy4)
    mesh_data.xiz .= Float64.(dummy4)
    read(f, dummy4)
    mesh_data.etax .= Float64.(dummy4)
    read(f, dummy4)
    mesh_data.etay .= Float64.(dummy4)
    read(f, dummy4)
    mesh_data.etaz .= Float64.(dummy4)
    read(f, dummy4)
    mesh_data.gammax .= Float64.(dummy4)
    read(f, dummy4)
    mesh_data.gammay .= Float64.(dummy4)
    read(f, dummy4)
    mesh_data.gammaz .= Float64.(dummy4)

    close(f)

    # separate crustal mesh layers for REGIONAL_MOHO_MESH 
    # 3-layer crust: 10(third layer), 11, 12(shallowest layer)
    if REGIONAL_MOHO_MESH
        num = 0
        for ispec = 1:nspec
            if mesh_data.idoubling[ispec] == IFLAG_CRUST
                id = num - fld(num, 3) * 3
                mesh_data.idoubling[ispec] = 10 * IFLAG_CRUST + id 
                num = num + 1
            end
        end
    end

    # separate mesh layers across 410-km
    # 40: above 410, 41: below 410
    for ispec = 1:nspec
        if mesh_data.idoubling[ispec] == IFLAG_670_220
            iglob = mesh_data.ibool[MIDX, MIDY, MIDZ, ispec]
            # element center coordinate
            xyz_center = mesh_data.xyz_glob[:,iglob]
            depth = (1.0 - sqrt(sum(mesh_data.xyz_glob[:,iglob].^2))) * R_EARTH_KM
            # this is dangerous due to 410 undulation
            # depth < 410 ? mesh_data.idoubling(ispec) = 10 * IFLAG_670_220 : mesh_data.idoubling(ispec) = 10 * IFLAG_670_220 + 1
            if depth < 410
                mesh_data.idoubling[ispec] = 10 * IFLAG_670_220
            else
                mesh_data.idoubling[ispec] = 10 * IFLAG_670_220 + 1
            end
        end
    end
    return mesh_data
end

function sem_io_open_file_for_read(basedir::String, iproc::Int64, tag::String)
    filename = "$(basedir)/proc$(lpad(string(iproc), 6, '0'))_reg1_$(tag).bin"
    f = FortranFile(filename)
    return f
end

function sem_io_open_file_for_write(basedir::String, iproc::Int64, tag::String)
    filename = "$(basedir)/proc$(lpad(string(iproc), 6, '0'))_reg1_$(tag).bin"
    f = FortranFile(filename, "w")
    return f
end

function sem_mesh_init!(mesh_data::sem_mesh_data, nspec::Int64, nglob::Int64)
    mesh_data.xyz_glob = zeros(Float64, 3, nglob)
    mesh_data.ibool = zeros(Int64, NGLLX, NGLLY, NGLLZ, nspec)
    mesh_data.idoubling = zeros(Int64, nspec)
    mesh_data.ispec_is_tiso = zeros(Bool, nspec)
    mesh_data.xix = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    mesh_data.xiy = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    mesh_data.xiz = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    mesh_data.etax = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    mesh_data.etay = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    mesh_data.etaz = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    mesh_data.gammax = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    mesh_data.gammay = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    mesh_data.gammaz = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
end

# * data 
function sem_io_read_gll_file_1!(basedir::String, iproc::Int64, model_name::String, model_gll::Array{Float64,4})
    nspec = size(model_gll)[4]
    f = sem_io_open_file_for_read(basedir, iproc, model_name)

    dummy = zeros(Float32, NGLLX, NGLLY, NGLLZ, nspec)
    read(f, dummy)
    close(f)
    model_gll .= Float64.(dummy)
end

function sem_io_read_gll_file_n!(basedir::String, iproc::Int64, model_names::Vector{String}, nmodel::Int64, model_gll::Array{Float64,5})
    dummy = similar(model_gll[1,:,:,:,:])
    for imodel = 1:nmodel
        sem_io_read_gll_file_1!(basedir, iproc, model_names[imodel], dummy)
        model_gll[imodel,:,:,:,:] = dummy
    end
end

function sem_io_write_gll_file_n(basedir::String, iproc::Int64, model_names::Vector{String}, nmodel::Int64, model_gll::Array{Float64,5})
    for imodel = 1:nmodel
        sem_io_write_gll_file_1(basedir, iproc, model_names[imodel], model_gll[imodel,:,:,:,:])
    end
end

function sem_io_write_gll_file_1(basedir::String, iproc::Int64, model_name::String, model_gll::Array{Float64,4})
    f = sem_io_open_file_for_write(basedir, iproc, model_name)
    dummy = Float32.(model_gll)
    write(f, dummy)
    close(f)
end