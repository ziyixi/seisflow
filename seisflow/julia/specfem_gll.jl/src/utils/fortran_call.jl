const file_path = Base.source_path()
const file_dir = splitdir(file_path)[1]

const lib_gll_library = joinpath(file_dir, "../fortran/lib/gll_library.so")
const lib_sem_mesh_mod = joinpath(file_dir, "../fortran/lib/sem_mesh_mod.so")

# gll_library
function zwgljd!(z::Vector{Float64}, w::Vector{Float64}, np::Int64, alpha::Float64, beta::Float64)
    # convert to references
    np_ref = Ref{Int32}(Int32(np))
    alpha_ref = Ref{Float64}(alpha)
    beta_ref = Ref{Float64}(beta)

    ccall((:zwgljd_, lib_gll_library), Cvoid, (Ptr{Float64}, Ptr{Float64}, Ref{Int32}, Ref{Float64}, Ref{Float64}), z, w, np_ref, alpha_ref, beta_ref)
end

# sem_mesh_mod
function anchor_point_index!(iax::Vector{Int32}, iay::Vector{Int32}, iaz::Vector{Int32})
    iax_32 = Int32.(iax)
    iay_32 = Int32.(iay)
    iaz_32 = Int32.(iaz)
    @assert length(iax) == 27 && length(iay) == 27 && length(iaz) == 27
    ccall((:anchor_point_index_, lib_sem_mesh_mod), Cvoid, (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}), iax_32, iay_32, iaz_32)
    iax .= Int32.(iax_32)
    iay .= Int32.(iay_32)
    iaz .= Int32.(iaz_32)
end

# sem_mesh_mod
function xyz2cube_bounded!(xyz_anchor::Array{Float64,2}, xyz::Array{Float64,1}, uvw::Array{Float64,1}, misloc::Float64, flag_inside::Bool)
    misloc_ref = Ref{Float64}(misloc)
    flag_inside_ref = Ref{Bool}(flag_inside)
    ccall((:xyz2cube_bounded_, lib_sem_mesh_mod), Cvoid, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}, Ref{Bool}), xyz_anchor, xyz, uvw, misloc_ref, flag_inside_ref)
    return misloc_ref.x, flag_inside_ref.x
end

# sem_mesh_mod
function lagrange_poly!(x::Float64, ngll::Int64, xgll::Vector{Float64}, lagrange::Vector{Float64})
    x_ref = Ref{Float64}(x)
    ngll_ref = Ref{Int32}(Int32(ngll))
    ccall((:lagrange_poly_, lib_sem_mesh_mod), Cvoid, (Ref{Float64}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}), x_ref, ngll_ref, xgll, lagrange)
end