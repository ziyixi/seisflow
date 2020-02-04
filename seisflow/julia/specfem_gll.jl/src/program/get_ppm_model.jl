import MPI

include("../types/types.jl")
include("../setting/constants.jl")
include("../utils/readfiles.jl")
include("../utils/kdtree.jl")
include("../utils/parse_commandline_ppm.jl")

function array_split_mpi(input_array,num_split,rank)
    div_result=div(length(input_array),num_split)
    mod_result=mod(length(input_array),num_split)
    startindex=nothing
    endindex=nothing
    if rank<=mod_result
        startindex=(div_result+1)*(rank-1)+1
        endindex=startindex+div_result
    else 
        startindex=(div_result+1)*mod_result+div_result*(rank-mod_result-1)+1
        endindex=startindex+div_result-1
    end
    return input_array[startindex:endindex]
end


function generate_profile_points(latnpts, lonnpts, vnpts, lon1, lat1, dep1, lon2, lat2, dep2, comm, latproc, lonproc)
    # MPI parameters
    rank = MPI.Comm_rank(comm)
    size_mpi = MPI.Comm_size(comm)

    # get ranges for the three directions
    rank_lat = rank % latproc + 1
    rank_lon = div(rank, latproc) + 1
    coor_lat=array_split_mpi(1:latnpts,latproc,rank_lat)
    coor_lon=array_split_mpi(1:lonnpts,lonproc,rank_lon)

    # init rθϕ_new
    ngll_new_this_rank = length(coor_lat) * length(coor_lon) * vnpts

    deplatlon_new_this_rank = zeros(Float64, 3, ngll_new_this_rank)
    # fill in deplatlon_new of this rank
    latnpts_this_rank = length(coor_lat)
    lonnpts_this_rank = length(coor_lon)
    for (vindex, dep) in enumerate(range(dep1, stop = dep2, length = vnpts))
        for (latindex, lat) in enumerate(range(lat1, stop = lat2, length = latnpts)[coor_lat])
            for (lonindex, lon) in enumerate(range(lon1, stop = lon2, length = lonnpts)[coor_lon])
                id = (vindex - 1) * latnpts_this_rank * lonnpts_this_rank + (latindex - 1) * lonnpts_this_rank + lonindex
                deplatlon_new_this_rank[:,id] = [dep,lat,lon]
            end
        end
    end

    xyz_new = zeros(Float64, 3, ngll_new_this_rank)
    # convert dep,lat,lon to x,y,z
    for id in 1:ngll_new_this_rank
        (dep, lat, lon) = deplatlon_new_this_rank[:,id]
        R_EARTH_KM = 6371.e0
        
        # normalized r
        r = (R_EARTH_KM - dep) / R_EARTH_KM
        θ = 90 - lat
        ϕ = lon

        # get x,y,z
        z = r * cosd(θ)
        h = r * sind(θ)
        x = h * cosd(ϕ)
        y = h * sind(ϕ)

        xyz_new[:,id] = [x,y,z]
    end

    return xyz_new
end


function run_interp(command_args::Dict{String,Any}, comm::MPI.Comm)
    # MPI
    rank = MPI.Comm_rank(comm)
    isroot = (rank == 0)

    # * init some variables
    nproc_old = command_args["nproc_old"]
    old_mesh_dir = command_args["old_mesh_dir"]
    old_model_dir = command_args["old_model_dir"]
    model_tags = command_args["model_tags"]
    output_file = command_args["output_file"]
    region = command_args["region"]
    npts = command_args["npts"]
    nproc = command_args["nproc"]
    (lon1, lat1, lon2, lat2, dep1, dep2) = parse.(Float64, split(region, "/"))
    (latnpts, lonnpts, vnpts) = parse.(Int64, split(npts, "/"))
    (latproc, lonproc) = parse.(Int64, split(nproc, "/"))

    # * init some variables
    mesh_old = sem_mesh_data()

    # * parse model tags 
    model_names = String.(split(model_tags, ","))
    nmodel = length(model_names)

    if isroot
        @info "# nmodel=$(nmodel)"
        @info "# model names=$(model_names)"
    end

    # * interpolate old mesh/model onto new mesh

    # typical element size at surface for old mesh
    # read from constants.jl for old mesh
    typical_size = deg2rad(max(ANGULAR_WIDTH_XI_IN_DEGREES_VAL / NEX_XI_VAL, ANGULAR_WIDTH_ETA_IN_DEGREES_VAL / NEX_ETA_VAL)) * R_UNIT_SPHERE
    max_search_dist = 10.0 * typical_size
    max_misloc = typical_size / 4.0

    # get xyz_new
    xyz_new = generate_profile_points(latnpts, lonnpts, vnpts, lon1, lat1, dep1, lon2, lat2, dep2, comm, latproc, lonproc)

    ngll_new_this_rank = size(xyz_new)[2]
    idoubling_new = zeros(Int64, ngll_new_this_rank)
    idoubling_new .= IFLAG_DUMMY

    # * initialize variables for interpolation
    location_1slice = Array{sem_mesh_location}(undef, ngll_new_this_rank)
    stat_final = zeros(Int64, ngll_new_this_rank)
    misloc_final = zeros(Float64, ngll_new_this_rank)

    stat_final .= -1
    misloc_final .= HUGEVAL

    model_interp_this_rank = zeros(nmodel, ngll_new_this_rank)
    # make model_interp_this_rank as 99999 to remove the poitns that hasn't been interpolated
    model_interp_this_rank .= 99999.0

    # * loop each slices of the old mesh
    flag = true

    # * loop for all points in xyz_new
    for iproc_old = 0:nproc_old - 1
        @info "[$rank]# iproc_old=$(iproc_old)"
        # read old mesh slice
        mesh_old = sem_mesh_read(old_mesh_dir, iproc_old)
        nspec_old = mesh_old.nspec
        # test if the new and old mesh slices are separated apart
        min_dist = HUGEVAL
        for ispec = 1:nspec_old
            iglob = mesh_old.ibool[MIDX,MIDY,MIDZ,ispec]
            old_x = mesh_old.xyz_glob[1,iglob]
            old_y = mesh_old.xyz_glob[2,iglob]       
            old_z = mesh_old.xyz_glob[3,iglob]
            dist_this_spec = sqrt(minimum(@. (xyz_new[1,:] - old_x)^2 + (xyz_new[2,:] - old_y)^2 + (xyz_new[3,:] - old_z)^2))
            min_dist = min(min_dist, dist_this_spec)
        end
        if min_dist > max_search_dist
            @info "[$rank]#$(iproc_old) slices too far away, skip"
            continue
        end

        # read old model
        model_gll_old = zeros(Float64, nmodel, NGLLX, NGLLY, NGLLZ, nspec_old)
        sem_io_read_gll_file_n!(old_model_dir, iproc_old, model_names, nmodel, model_gll_old)

        # locate points in this mesh slice
        nnearest = 10
        location_1slice = sem_mesh_locate_kdtree2!(mesh_old, ngll_new_this_rank, xyz_new, idoubling_new, nnearest, max_search_dist, max_misloc, iproc_old)

        for igll = 1:ngll_new_this_rank
            if (stat_final[igll] == 1 && location_1slice[igll].stat == 1)
                @info "[$rank]# multi-located, $(xyz_new[:,igll])"
                continue
            end
            # for point located inside one element in the first time or closer to one element than located before
            if location_1slice[igll].stat == 1 || (location_1slice[igll].stat == 0 && location_1slice[igll].misloc < misloc_final[igll])
                for imodel = 1:nmodel
                    model_interp_this_rank[imodel,igll] = sum(location_1slice[igll].lagrange .* model_gll_old[imodel,:,:,:,location_1slice[igll].eid])
                end
                stat_final[igll] = location_1slice[igll].stat
                misloc_final[igll] = location_1slice[igll].misloc
            end
        end
    end

    # gather all ngll_new_this_rank to ngll_new


    # * write out gll files for this new mesh slice
    # get latnpts_this_rank and lonnpts_this_rank
    rank = MPI.Comm_rank(comm)
    size_mpi = MPI.Comm_size(comm)

    rank_lat = rank % latproc + 1
    rank_lon = div(rank, latproc) + 1
    coor_lat = np[:array_split](1:latnpts, latproc)[rank_lat]
    coor_lon = np[:array_split](1:lonnpts, lonproc)[rank_lon]

    latnpts_this_rank = length(coor_lat)
    lonnpts_this_rank = length(coor_lon)

    # output files
    if isroot
        run(`rm -rf $output_file`)
        run(`mkdir -p $output_file`)
    end
    MPI.Barrier(comm)

    open(output_file * "/$rank", "w") do io
        for (vindex, dep) in enumerate(range(dep1, stop = dep2, length = vnpts))
            for (latindex, lat) in enumerate(range(lat1, stop = lat2, length = latnpts)[coor_lat])
                for (lonindex, lon) in enumerate(range(lon1, stop = lon2, length = lonnpts)[coor_lon])
                    id = (vindex - 1) * latnpts_this_rank * lonnpts_this_rank + (latindex - 1) * lonnpts_this_rank + lonindex
                    thesize = size(model_interp_this_rank)[1]
                    if thesize == 1
                        write(io, "$lon $lat $dep $(model_interp_this_rank[1,id]) \n")
                    else
                        out_string = "$lon $lat $dep "
                        for item in model_interp_this_rank[:,id]
                            out_string *= "$item "
                        end
                        out_string *= "\n"
                        write(io, out_string)
                    end
                end
            end
        end
    end

end


function main()
    MPI.Init()
    comm = MPI.COMM_WORLD

    command_args = parse_commandline()  
    run_interp(command_args, comm)  

    MPI.Finalize()
end

main()