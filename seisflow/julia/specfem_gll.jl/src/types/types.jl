mutable struct sem_mesh_data
    nspec::Int64
    nglob::Int64
    xyz_glob::Array{Float64,2}
    ibool::Array{Int64,4}
    idoubling::Array{Int64,1}
    ispec_is_tiso::Array{Bool,1}

    # xix::Array{Float64,4}
    # xiy::Array{Float64,4}
    # xiz::Array{Float64,4}
    # etax::Array{Float64,4}
    # etay::Array{Float64,4}
    # etaz::Array{Float64,4}
    # gammax::Array{Float64,4}
    # gammay::Array{Float64,4}
    # gammaz::Array{Float64,4}

    sem_mesh_data() = new()
end  

mutable struct sem_mesh_location
    stat::Int64
    eid::Int64

    uvw::Vector{Float64}
    misloc::Float64
    lagrange::Array{Float64,3}

    sem_mesh_location() = new()
end

mutable struct sem_boundary_disc
    NSPEC2D_MOHO::Int64
    NSPEC2D_400::Int64
    NSPEC2D_670::Int64
    ibelm_moho_top::Array{Int64,1}
    ibelm_moho_bot::Array{Int64,1}
    ibelm_400_top::Array{Int64,1}
    ibelm_400_bot::Array{Int64,1}
    ibelm_670_top::Array{Int64,1}
    ibelm_670_bot::Array{Int64,1}

    sem_boundary_disc() = new()
end