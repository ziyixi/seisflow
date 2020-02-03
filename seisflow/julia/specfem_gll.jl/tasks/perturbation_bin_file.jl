include("../src/utils/readfiles.jl")
include("../src/setting/constants.jl")

"""
generate perturbation file according to the reference bin file
"""
function generate_perturbation(target_basedir::String, reference_basedir::String, output_basedir::String, nproc::Int64, nspec::Int64)
    tags = ["vph","vpv","vsh","vsv","eta","qmu","rho"]
    model_gll_reference = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    model_gll_target = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    model_gll_output = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    for iproc in 0:nproc - 1
        if(iproc % 10 == 0)
            @info "step" iproc
        end
        for tag in tags
            sem_io_read_gll_file_1!(reference_basedir, iproc, tag, model_gll_reference)
            sem_io_read_gll_file_1!(target_basedir, iproc, tag, model_gll_target)
            
            model_gll_output = (model_gll_target .- model_gll_reference) ./ model_gll_reference
            @assert count(isnan.(model_gll_output)) == 0
            sem_io_write_gll_file_1(output_basedir, iproc, tag, model_gll_output)
        end
    end
end

"""
generate the real vaue bin file based on reference model and perturbation bin file
"""
function generate_real(target_basedir::String, reference_basedir::String, output_basedir::String, nproc::Int64, nspec::Int64)
    tags = ["vph","vpv","vsh","vsv","eta","qmu","rho"]
    model_gll_reference = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    model_gll_target = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    model_gll_output = zeros(Float64, NGLLX, NGLLY, NGLLZ, nspec)
    for iproc in 0:nproc - 1
        if(iproc % 10 == 0)
            @info "step" iproc
        end
        for tag in tags
            sem_io_read_gll_file_1!(reference_basedir, iproc, tag, model_gll_reference)
            sem_io_read_gll_file_1!(target_basedir, iproc, tag, model_gll_target)
            
            model_gll_output = model_gll_reference .* (1 .+ model_gll_target)
            @assert count(isnan.(model_gll_output)) == 0
            sem_io_write_gll_file_1(output_basedir, iproc, tag, model_gll_output)
        end
    end
end