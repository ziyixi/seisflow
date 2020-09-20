include("../src/utils/readfiles.jl")
include("../src/setting/constants.jl")
using ProgressMeter
using Base.Threads

"""
generate perturbation file according to the reference bin file
"""
function generate_perturbation(target_basedir::String, reference_basedir::String, output_basedir::String, tags::String, nproc::Int64, nspec::Int64)
    # p = Progress(nproc)
    @threads for iproc in 0:nproc - 1
        kernel_generate_perturbation(target_basedir, reference_basedir, output_basedir, tags, iproc, nspec)
        # next!(p)
    end
end

function kernel_generate_perturbation(target_basedir::String, reference_basedir::String, output_basedir::String, tags::String, iproc::Int64, nspec::Int64)
    tags_splitted = split(tags, ",")
    model_gll_reference = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    model_gll_target = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    model_gll_output = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    for tag in tags_splitted
        # convert tag to String
        tag = String(tag)
        sem_io_read_gll_file_1!(reference_basedir, iproc, tag, model_gll_reference)
        sem_io_read_gll_file_1!(target_basedir, iproc, tag, model_gll_target)
            
        model_gll_output = (model_gll_target .- model_gll_reference) ./ model_gll_reference
        sem_io_write_gll_file_1(output_basedir, iproc, tag, model_gll_output)
    end
end

"""
generate the real vaue bin file based on reference model and perturbation bin file
"""
function generate_real(target_basedir::String, reference_basedir::String, output_basedir::String, tags::String, nproc::Int64, nspec::Int64)
    tags_splitted = split(tags, ",")
    p = Progress(nproc)
    @threads for iproc in 0:nproc - 1
        model_gll_reference = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
        model_gll_target = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
        model_gll_output = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
        for tag in tags_splitted
            # convert tag to String
            tag = String(tag)
            sem_io_read_gll_file_1!(reference_basedir, iproc, tag, model_gll_reference)
            sem_io_read_gll_file_1!(target_basedir, iproc, tag, model_gll_target)
            
            model_gll_output = model_gll_reference .* (1 .+ model_gll_target)
            sem_io_write_gll_file_1(output_basedir, iproc, tag, model_gll_output)
        end
        next!(p)
    end
end

"""
calculate the difference between two sets of gll files.
"""
function generate_difference(target_basedir::String, reference_basedir::String, output_basedir::String, tags::String, nproc::Int64, nspec::Int64)
    tags_splitted = split(tags, ",")
    p = Progress(nproc)
    @threads for iproc in 0:nproc - 1
        model_gll_reference = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
        model_gll_target = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
        model_gll_output = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
        for tag in tags_splitted
            # convert tag to String
            tag = String(tag)
            sem_io_read_gll_file_1!(reference_basedir, iproc, tag, model_gll_reference)
            sem_io_read_gll_file_1!(target_basedir, iproc, tag, model_gll_target)
            
            model_gll_output = model_gll_target .- model_gll_reference
            sem_io_write_gll_file_1(output_basedir, iproc, tag, model_gll_output)
        end
        next!(p)
    end
end

"""
calculate the model with the minus sign
"""
function minus_sign(target_basedir::String, output_basedir::String, tags::String, nproc::Int64, nspec::Int64)
    tags_splitted = split(tags, ",")
    p = Progress(nproc)
    @threads for iproc in 0:nproc - 1
        model_gll_target = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
        model_gll_output = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
        for tag in tags_splitted
            # convert tag to String
            tag = String(tag)
            sem_io_read_gll_file_1!(target_basedir, iproc, tag, model_gll_target)
            
            model_gll_output = -model_gll_target
            sem_io_write_gll_file_1(output_basedir, iproc, tag, model_gll_output)
        end
        next!(p)
    end
end

"""
get new model from the perturbation. 
"""
function generate_new(old_basedir::String, per_basedir::String, output_basedir::String, tags::String, nproc::Int64, nspec::Int64)
    p = Progress(nproc)
    @threads for iproc in 0:nproc - 1
        kernel_generate_new(old_basedir, per_basedir, output_basedir, tags, iproc, nspec)
        next!(p)
    end
end

function kernel_generate_new(old_basedir::String, per_basedir::String, output_basedir::String, tags::String, iproc::Int64, nspec::Int64)
    tags_splitted = split(tags, ",")
    model_gll_per = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    model_gll_old = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    model_gll_output = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    for tag in tags_splitted
        try
            # convert tag to String
            tag = String(tag)
            tag_old = tag * "_new"
            sem_io_read_gll_file_1!(old_basedir, iproc, tag_old, model_gll_old)
            # modify tag here
            tag_per = "d" * tag * tag
            sem_io_read_gll_file_1!(per_basedir, iproc, tag_per, model_gll_per)
                
            model_gll_output = model_gll_old .* (exp.(model_gll_per))
            sem_io_write_gll_file_1(output_basedir, iproc, tag_old, model_gll_output)
        catch e
        end
    end
end

"""
Compare two models.
"""
function compare_model(target_basedir::String, reference_basedir::String, tag::String,  nspec::Int64, iproc::Int64)
    model_gll_target = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    model_gll_reference = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    sem_io_read_gll_file_1!(target_basedir, iproc, tag, model_gll_target)
    sem_io_read_gll_file_1!(reference_basedir, iproc, tag, model_gll_reference)
    model_diff = model_gll_target .- model_gll_reference
    @info "max difference" maximum(model_diff)
end

"""
convert vpv,vph,vsv,vsh to bulk_cv,bulk_ch,vsv,vsh.
"""
function vpvs2bulk_c!(target_basedir::String, output_basedir::String, nspec::Int64, iproc::Int64)
    model_gll_vpv = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    model_gll_vph = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    model_gll_vsv = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    model_gll_vsh = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    # * read in models
    sem_io_read_gll_file_1!(target_basedir, iproc, "vpv", model_gll_vpv)
    sem_io_read_gll_file_1!(target_basedir, iproc, "vph", model_gll_vph)
    sem_io_read_gll_file_1!(target_basedir, iproc, "vsv", model_gll_vsv)
    sem_io_read_gll_file_1!(target_basedir, iproc, "vsh", model_gll_vsh)
    # * now we override vpv,vph with the bulk_c, calculate based on:
    # * ϕ^2=α^2-4/3*β^2
    model_gll_vpv = @. sqrt(model_gll_vpv^2 - 4 / 3 * model_gll_vsv^2)
    model_gll_vph = @. sqrt(model_gll_vph^2 - 4 / 3 * model_gll_vsh^2)
    # * write out files
    sem_io_write_gll_file_1(output_basedir, iproc, "bulk_cv", model_gll_vpv)
    sem_io_write_gll_file_1(output_basedir, iproc, "bulk_ch", model_gll_vph)
    sem_io_write_gll_file_1(output_basedir, iproc, "vsv", model_gll_vsv)
    sem_io_write_gll_file_1(output_basedir, iproc, "vsh", model_gll_vsh)
end

"""
convert bulk_cv,bulk_ch,vsv,vsh to vpv,vph,vsv,vsh.
"""
function bulk_c2vpvs!(target_basedir::String, output_basedir::String, nspec::Int64, iproc::Int64)
    model_gll_vpv = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    model_gll_vph = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    model_gll_vsv = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    model_gll_vsh = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    # * read in models
    sem_io_read_gll_file_1!(target_basedir, iproc, "bulk_cv", model_gll_vpv)
    sem_io_read_gll_file_1!(target_basedir, iproc, "bulk_ch", model_gll_vph)
    sem_io_read_gll_file_1!(target_basedir, iproc, "vsv", model_gll_vsv)
    sem_io_read_gll_file_1!(target_basedir, iproc, "vsh", model_gll_vsh)
    # * now we override bilk_c with vpv and vph, calculate based on:
    # * ϕ^2=α^2-4/3*β^2 so α^2=ϕ^2+4/3*β^2
    model_gll_vpv = @. sqrt(model_gll_vpv^2 + 4 / 3 * model_gll_vsv^2)
    model_gll_vph = @. sqrt(model_gll_vph^2 + 4 / 3 * model_gll_vsh^2)
    # * write out files
    sem_io_write_gll_file_1(output_basedir, iproc, "vpv", model_gll_vpv)
    sem_io_write_gll_file_1(output_basedir, iproc, "vph", model_gll_vph)
    sem_io_write_gll_file_1(output_basedir, iproc, "vsv", model_gll_vsv)
    sem_io_write_gll_file_1(output_basedir, iproc, "vsh", model_gll_vsh)
end