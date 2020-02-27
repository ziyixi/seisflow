include("../src/utils/readfiles.jl")
include("../src/setting/constants.jl")
using ProgressMeter
using Base.Threads

"""
generate perturbation file according to the reference bin file
"""
function generate_perturbation(target_basedir::String, reference_basedir::String, output_basedir::String, tags::String, nproc::Int64, nspec::Int64)
    p = Progress(nproc)
    @threads for iproc in 0:nproc - 1
        kernel_generate_perturbation(target_basedir, reference_basedir, output_basedir, tags, iproc, nspec)
        next!(p)
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
        try
            kernel_generate_new(old_basedir, per_basedir, output_basedir, tags, iproc, nspec)
        catch e
        end
        next!(p)
    end
end

function kernel_generate_new(old_basedir::String, per_basedir::String, output_basedir::String, tags::String, iproc::Int64, nspec::Int64)
    tags_splitted = split(tags, ",")
    model_gll_per = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    model_gll_old = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    model_gll_output = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    for tag in tags_splitted
        # convert tag to String
        tag = String(tag)
        tag_old=tag*"_new"
        sem_io_read_gll_file_1!(old_basedir, iproc, tag_old, model_gll_old)
        # modify tag here
        tag_per="d"*tag*tag
        sem_io_read_gll_file_1!(per_basedir, iproc, tag_per, model_gll_per)
            
        model_gll_output = model_gll_old.*(exp.(model_gll_per))
        sem_io_write_gll_file_1(output_basedir, iproc, tag_old, model_gll_output)
    end
end

"""
Compare two models.
"""
function compare_model(target_basedir::String, reference_basedir::String,tag::String,  nspec::Int64,iproc::Int64)
    model_gll_target = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    model_gll_reference = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    sem_io_read_gll_file_1!(target_basedir, iproc, tag, model_gll_target)
    sem_io_read_gll_file_1!(reference_basedir, iproc, tag, model_gll_reference)
    model_diff=model_gll_target.-model_gll_reference
    @info "max difference" maximum(model_diff)
end