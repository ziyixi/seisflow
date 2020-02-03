using NearestNeighbors
include("../types/types.jl")
include("../setting/constants.jl")
include("./fortran_call.jl")

function sem_mesh_locate_kdtree2!(mesh_data::sem_mesh_data, npoint::Int64, xyz::Array{Float64,2}, idoubling::Vector{Int64}, nnearest::Int64, max_search_dist::Float64, max_misloc::Float64, iproc_old::Int64)
    # * init result
    location_result = Vector{sem_mesh_location}(undef, npoint)
    for i in 1:npoint
        location_result[i] = sem_mesh_location()
    end

    # * GLL colocation points and lagrange interpolation weights
    xigll = zeros(Float64, NGLLX)
    wxgll = zeros(Float64, NGLLX)
    yigll = zeros(Float64, NGLLY)
    wygll = zeros(Float64, NGLLY)
    zigll = zeros(Float64, NGLLZ)
    wzgll = zeros(Float64, NGLLZ)
    hlagx = zeros(Float64, NGLLX)
    hlagy = zeros(Float64, NGLLY)
    hlagz = zeros(Float64, NGLLZ)

    # * mesh dimension variable
    nspec = mesh_data.nspec

    # * coordinates of GLL points
    zwgljd!(xigll, wxgll, NGLLX, GAUSSALPHA, GAUSSBETA)
    zwgljd!(yigll, wygll, NGLLY, GAUSSALPHA, GAUSSBETA)
    zwgljd!(zigll, wzgll, NGLLZ, GAUSSALPHA, GAUSSBETA)

    # * anchor points of mesh elements 

    # get index of anchor points in the GLL element
    iax = zeros(Int32, NGNOD)
    iay = zeros(Int32, NGNOD)
    iaz = zeros(Int32, NGNOD)
    anchor_point_index!(iax, iay, iaz)

    # get anchor and center points of each GLL element 
    xyz_elem = zeros(Float64, 3, nspec)
    xyz_anchor = zeros(Float64, 3, NGNOD, nspec)

    for ispec = 1:nspec
        for ia = 1:NGNOD
            iglob = mesh_data.ibool[iax[ia], iay[ia], iaz[ia], ispec]
            xyz_anchor[:,ia,ispec] = mesh_data.xyz_glob[:,iglob]
        end
        # the last anchor point is the element center
        xyz_elem[:,ispec] = xyz_anchor[:,NGNOD,ispec]
    end

    # * build kdtree
    kdtree = KDTree(xyz_elem)
    # * locate each point
    # initialize location results
    for i in 1:npoint
        location_result[i].stat = -1
        location_result[i].eid = -1
        location_result[i].misloc = HUGEVAL
        location_result[i].lagrange = zeros(Float64, NGLLX, NGLLY, NGLLZ)
        location_result[i].uvw = zeros(Float64, 3)
    end

    # loop points
    for ipoint = 1:npoint
        xyz1 = xyz[:,ipoint]
       # * get the n nearest elements in the mesh
        idxs, _ = knn(kdtree, xyz1, nnearest)

       # * test each neighbour elements to see if target point is located inside
        for inn = 1:nnearest
            ispec = idxs[inn]
            # skip the element a certain distance away
            dist = sqrt(sum((xyz_elem[:,ispec] .- xyz1).^2))
            # if (ipoint == 492835) 
            #     @info "@@@@"
            #     @info dist, max_search_dist, inn
            # end
            if dist > max_search_dist
                continue
            end

            # only use element with the same layer ID (i.e. idoubling)
            if (idoubling[ipoint] != IFLAG_DUMMY && mesh_data.idoubling[ispec] != idoubling[ipoint])
                continue
            end

            # locate point to this element
            uvw1 = similar(xyz1)
            misloc1 = HUGEVAL
            flag_inside = false
            misloc1, flag_inside = xyz2cube_bounded!(xyz_anchor[:,:,ispec], xyz1, uvw1, misloc1, flag_inside)
            # if (ipoint == 492835) 
            #     @info "@@@@"
            #     @info flag_inside, misloc1, uvw1, ispec
            # end
            if flag_inside == true
                location_result[ipoint].stat = 1
                location_result[ipoint].eid = ispec
                location_result[ipoint].misloc = misloc1
                location_result[ipoint].uvw = uvw1
                break
            else
                if (misloc1 < max_misloc) && (misloc1 < location_result[ipoint].misloc)
                    location_result[ipoint].stat = 0
                    location_result[ipoint].eid = ispec
                    location_result[ipoint].misloc = misloc1
                    location_result[ipoint].uvw = uvw1
                end
            end
        end

        # * set interpolation weights on GLL points if located
        if location_result[ipoint].stat != -1
            lagrange_poly!(location_result[ipoint].uvw[1], NGLLX, xigll, hlagx)
            lagrange_poly!(location_result[ipoint].uvw[2], NGLLY, yigll, hlagy)
            lagrange_poly!(location_result[ipoint].uvw[3], NGLLZ, zigll, hlagz)

            for igllz = 1:NGLLZ
                for iglly = 1:NGLLY
                    for igllx = 1:NGLLX
                        location_result[ipoint].lagrange[igllx,iglly,igllz] = hlagx[igllx] * hlagy[iglly] * hlagz[igllz]
                    end
                end
            end
        end
        #to debug
        # if (ipoint == 492835) && (location_result[ipoint].stat != -1)
        #     @info "@@@@"
        #     @info iproc_old
        #     @info xyz1
        #     @info idxs
        #     @info hlagx
        #     @info hlagy
        #     @info hlagz
        #     @info location_result[ipoint].uvw
        #     @info "@@@@"
        # end
    end  
    return location_result
end