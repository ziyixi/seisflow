import MPI

include("../types/types.jl")
include("../setting/constants.jl")
include("../utils/readfiles.jl")
include("../utils/kdtree.jl")
include("../utils/parse_commandline.jl")

# ! this program use Float64, which may cause the out of memory problem when In values_from_mesh.h 
# ! approximate total number of points in entire mesh > 151837623.000000.
function run_interp(myrank::Int64, nrank::Int64, command_args::Dict{String,Any})
    # * init some variables
    isroot = (myrank == 0)
    nproc_old = command_args["nproc_old"]
    old_mesh_dir = command_args["old_mesh_dir"]
    old_model_dir = command_args["old_model_dir"]
    nproc_new = command_args["nproc_new"]
    new_mesh_dir = command_args["new_mesh_dir"]
    new_model_dir = command_args["new_model_dir"]
    model_tags = command_args["model_tags"]
    output_dir = command_args["output_dir"]

    # * init some variables
    mesh_new = sem_mesh_data()
    mesh_old = sem_mesh_data()

    # * parse model tags 
    model_names = String.(split(model_tags, ","))
    nmodel = length(model_names)
    if isroot
        @info "[$(myrank)]# nmodel=$(nmodel)"
        @info "[$(myrank)]# model names=$(model_names)"
    end


    # * interpolate old mesh/model onto new mesh

    # typical element size at surface for old mesh
    # read from constants.jl for old mesh
    typical_size = deg2rad(max(ANGULAR_WIDTH_XI_IN_DEGREES_VAL / NEX_XI_VAL, ANGULAR_WIDTH_ETA_IN_DEGREES_VAL / NEX_ETA_VAL)) * R_UNIT_SPHERE
    max_search_dist = 10.0 * typical_size
    max_misloc = typical_size / 4.0

    # * loop each new mesh slice (where MPI parallels)
    for iproc_new in myrank:nrank:(nproc_new - 1)
        # if iproc_new != 386
        #     continue
        # end
        @info "[$(myrank)]# iproc_new=$(iproc_new)"

        # * read new mesh slice
        mesh_new = sem_mesh_read(new_mesh_dir, iproc_new)
        nspec_new = mesh_new.nspec
        ngll_new = NGLLX * NGLLY * NGLLZ * nspec_new

        xyz_new = zeros(Float64, 3, ngll_new)
        idoubling_new = zeros(Int64, ngll_new)
        xyz_center_new = zeros(Float64, 3, nspec_new)

        for ispec = 1:nspec_new
            iglob = mesh_new.ibool[MIDX,MIDY,MIDZ,ispec]
            xyz_center_new[:,ispec] = mesh_new.xyz_glob[:,iglob]
            for igllz = 1:NGLLZ
                for iglly = 1:NGLLY
                    for igllx = 1:NGLLX
                        igll = igllx + NGLLX * ( (iglly - 1) + NGLLY * ( (igllz - 1) + NGLLZ * ( (ispec - 1))))
                        iglob = mesh_new.ibool[igllx,iglly,igllz,ispec]
                        xyz_new[:,igll] = mesh_new.xyz_glob[:,iglob]
                        idoubling_new[igll] = mesh_new.idoubling[ispec]
                    end
                end
            end
        end

        # * initialize variables for interpolation
        location_1slice = Array{sem_mesh_location}(undef, ngll_new)
        stat_final = zeros(Int64, ngll_new)
        misloc_final = zeros(Float64, ngll_new)

        stat_final .= -1
        misloc_final .= HUGEVAL

        model_interp = zeros(nmodel, ngll_new)

        # * read in the new model as the background model
        model_gll_new = zeros(Float64, nmodel, NGLLX, NGLLY, NGLLZ, nspec_new)
        sem_io_read_gll_file_n!(new_model_dir, iproc_new, model_names, nmodel, model_gll_new)

        # * loop each slices of the old mesh
        flag = true
        for iproc_old = 0:nproc_old - 1
            @info "[$(myrank)]# iproc_old=$(iproc_old)"
            # read old mesh slice
            mesh_old = sem_mesh_read(old_mesh_dir, iproc_old)
            nspec_old = mesh_old.nspec
            # test if the new and old mesh slices are separated apart
            min_dist = HUGEVAL
            for ispec = 1:nspec_old
                iglob = mesh_old.ibool[MIDX,MIDY,MIDZ,ispec]
                min_dist = min(min_dist, sqrt(minimum((xyz_center_new[1,:] .- mesh_old.xyz_glob[1,iglob]).^2 + (xyz_center_new[2,:] .- mesh_old.xyz_glob[2,iglob]).^2 + (xyz_center_new[3,:] .- mesh_old.xyz_glob[3,iglob]).^2)))
            end
            if min_dist > max_search_dist
                @info "[$(myrank)]# $(iproc_old)/$(iproc_new) slices too far away, skip"
                continue
            end

            # read old model
            model_gll_old = zeros(Float64, nmodel, NGLLX, NGLLY, NGLLZ, nspec_old)
            sem_io_read_gll_file_n!(old_model_dir, iproc_old, model_names, nmodel, model_gll_old)

            # locate points in this mesh slice
            nnearest = 10
            location_1slice = sem_mesh_locate_kdtree2!(mesh_old, ngll_new, xyz_new, idoubling_new, nnearest, max_search_dist, max_misloc, iproc_old)
            for igll = 1:ngll_new
                if (stat_final[igll] == 1 && location_1slice[igll].stat == 1)
                    @info "[$(myrank)]# multi-located, $(xyz_new[:,igll])"
                    continue
                end
                # for point located inside one element in the first time or closer to one element than located before
                if location_1slice[igll].stat == 1 || (location_1slice[igll].stat == 0 && location_1slice[igll].misloc < misloc_final[igll])
                    for imodel = 1:nmodel
                        model_interp[imodel,igll] = sum(location_1slice[igll].lagrange .* model_gll_old[imodel,:,:,:,location_1slice[igll].eid])
                    end
                    stat_final[igll] = location_1slice[igll].stat
                    misloc_final[igll] = location_1slice[igll].misloc
                end
                # if location_1slice[igll].stat == 1 || (location_1slice[igll].stat == 0 && location_1slice[igll].misloc < misloc_final[igll])
                    # flag = false
                # end
            end
        end
        # * write out gll files for this new mesh slice

        # reshape model_interp to model_gll
        for ispec = 1:nspec_new
            for igllz = 1:NGLLZ
                for iglly = 1:NGLLY
                    for igllx = 1:NGLLX
                        igll = igllx + NGLLX * ( (iglly - 1) + NGLLY * ( (igllz - 1) + NGLLZ * ( (ispec - 1))))
                        if stat_final[igll] != -1
                            model_gll_new[:,igllx,iglly,igllz,ispec] = model_interp[:,igll]
                        end
                    end
                end
            end
        end
        sem_io_write_gll_file_n(output_dir, iproc_new, model_names, nmodel, model_gll_new)
    end
end
    




function main()
    MPI.Init()

    comm = MPI.COMM_WORLD
    myrank = MPI.Comm_rank(comm)
    nrank = MPI.Comm_size(comm)
    command_args = parse_commandline()
    
    run_interp(myrank, nrank, command_args)
    
    MPI.Barrier(comm)
    MPI.Finalize()
end

main()