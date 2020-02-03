using Geodesics

include("../types/types.jl")
include("../setting/constants.jl")
include("../utils/readfiles.jl")
include("../utils/kdtree.jl")
include("../utils/parse_commandline_profile.jl")


function generate_profile_points(hnpts, vnpts, lon1, lat1, dep1, lon2, lat2, dep2)
    # if cut vertically
    is_cut_vertically = true
    if lon1 == lon2
        is_cut_vertically = false
    end

    # init rθϕ_new
    ngll_new = hnpts * vnpts
    xyz_new = zeros(Float64, 3, ngll_new)
    deplatlon_new = zeros(Float64, 3, ngll_new)

    # fill in deplatlon_new
    if is_cut_vertically
        for (vindex, dep) in enumerate(range(dep1, stop = dep2, length = vnpts))
            for (hindex, (lon, lat)) in enumerate(zip(range(lon1, stop = lon2, length = hnpts), range(lat1, stop = lat2, length = hnpts)))
                id = (vindex - 1) * hnpts + hindex
                deplatlon_new[:,id] = [dep,lat,lon]
            end
        end
    else
        lonnpts = hnpts
        latnpts = vnpts 
        for (lonindex, lon) in enumerate(range(lon1, stop = lon2, length = lonnpts))
            for (latindex, lat) in enumerate(range(lat1, stop = lat2, length = latnpts))
                id = (lonindex - 1) * latnpts + latindex
                deplatlon_new[:,id] = [dep1,lat,lon]
            end
        end
    end

    # convert dep,lat,lon to x,y,z
    for id in 1:ngll_new
        (dep, lat, lon) = deplatlon_new[:,id]
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


function run_interp(command_args::Dict{String,Any})
    # * init some variables
    nproc_old = command_args["nproc_old"]
    old_mesh_dir = command_args["old_mesh_dir"]
    old_model_dir = command_args["old_model_dir"]
    model_tags = command_args["model_tags"]
    output_file = command_args["output_file"]
    region = command_args["region"]
    npts = command_args["npts"]
    (lon1, lat1, lon2, lat2, dep1, dep2) = parse.(Float64, split(region, "/"))
    (hnpts, vnpts) = parse.(Int64, split(npts, "/"))

    # * init some variables
    mesh_old = sem_mesh_data()

    # * parse model tags 
    model_names = String.(split(model_tags, ","))
    nmodel = length(model_names)

    @info "# nmodel=$(nmodel)"
    @info "# model names=$(model_names)"


    # * interpolate old mesh/model onto new mesh

    # typical element size at surface for old mesh
    # read from constants.jl for old mesh
    typical_size = deg2rad(max(ANGULAR_WIDTH_XI_IN_DEGREES_VAL / NEX_XI_VAL, ANGULAR_WIDTH_ETA_IN_DEGREES_VAL / NEX_ETA_VAL)) * R_UNIT_SPHERE
    max_search_dist = 10.0 * typical_size
    max_misloc = typical_size / 4.0

    # get xyz_new
    xyz_new = generate_profile_points(hnpts, vnpts, lon1, lat1, dep1, lon2, lat2, dep2)

    ngll_new = size(xyz_new)[2]
    idoubling_new = zeros(Int64, ngll_new)
    idoubling_new .= IFLAG_DUMMY

    # * initialize variables for interpolation
    location_1slice = Array{sem_mesh_location}(undef, ngll_new)
    stat_final = zeros(Int64, ngll_new)
    misloc_final = zeros(Float64, ngll_new)

    stat_final .= -1
    misloc_final .= HUGEVAL

    model_interp = zeros(nmodel, ngll_new)

    # * loop each slices of the old mesh
    flag = true

    # * loop for all points in xyz_new
    for iproc_old = 0:nproc_old - 1
        @info "# iproc_old=$(iproc_old)"
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
            @info "#$(iproc_old) slices too far away, skip"
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
                @info "# multi-located, $(xyz_new[:,igll])"
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
        end
    end

    # * write out gll files for this new mesh slice
    is_cut_vertically = true
    if lon1 == lon2
        is_cut_vertically = false
    end

    open(output_file, "w") do io
        if is_cut_vertically
            a, f = Geodesics.EARTH_R_MAJOR_WGS84, Geodesics.F_WGS84;
            for (vindex, dep) in enumerate(range(dep1, stop = dep2, length = vnpts))
                for (hindex, (lon, lat)) in enumerate(zip(range(lon1, stop = lon2, length = hnpts), range(lat1, stop = lat2, length = hnpts)))
                    id = (vindex - 1) * hnpts + hindex
                    dist, az, baz = Geodesics.inverse(deg2rad.((lon1, lat1, lon, lat))..., a, f)
                    if size(model_interp)[1] == 1
                        write(io, "$lon $lat $dep $(model_interp[1,id]) $(dist / 1000.0) \n")
                    else
                        write(io, "$lon $lat $dep $(model_interp[:,id]) $(dist / 1000.0) \n")
                    end
                end
            end
        else
            lonnpts = hnpts
            latnpts = vnpts 
            for (lonindex, lon) in enumerate(range(lon1, stop = lon2, length = lonnpts))
                for (latindex, lat) in enumerate(range(lat1, stop = lat2, length = latnpts))
                    id = (lonindex - 1) * latnpts + latindex
                    if size(model_interp)[1] == 1
                        write(io, "$lon $lat $dep1 $(model_interp[1,id]) \n")
                    else
                        write(io, "$lon $lat $dep1 $(model_interp[:,id]) \n")
                    end
                end
            end
        end
    end
end


function main()
    command_args = parse_commandline()  
    run_interp(command_args)  
end

main()