files = ["gll_library.f90","sem_mesh_mod.f90"]
FF = "gfortran"

file_path = Base.source_path()
file_dir = splitdir(file_path)[1]

for file in files
    fname = join(split(file, ".")[1:end - 1], ".")
    command = `$(FF) -fPIC -shared $(file_dir)/src/$(fname).f90 -o $(file_dir)/lib/$(fname).so`
    run(command)
end