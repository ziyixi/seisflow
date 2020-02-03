from slurmpy import Slurm

# ! change setting first!!!
# per_s362ani_good + per_tao ->  per_s362ani_good_tao
nproc_old = 336
old_mesh_dir = "/scratch/05880/tg851791/work/generate_hybrid_v703/gll_work/control_file/tao"
old_model_dir = "/scratch/05880/tg851791/work/generate_hybrid_v703/gll_work/perturbation/per_tao"
nproc_new = 324
new_mesh_dir = "/scratch/05880/tg851791/work/generate_small_v703/specfem/s362ani_good/DATABASES_MPI"
new_model_dir = "/scratch/05880/tg851791/work/generate_small_v703/perturbation/per_s362ani_good"
model_tags = ",".join(["vph", "vpv", "vsh", "vsv", "eta", "qmu", "rho"])
output_dir = "/scratch/05880/tg851791/work/generate_small_v703/perturbation/per_s362ani_good_tao"

command1 = f"ibrun julia ../specfem_gll.jl/src/program/xsem_interp_mesh2.jl --nproc_old {nproc_old} --old_mesh_dir {old_mesh_dir} --old_model_dir {old_model_dir} --nproc_new {nproc_new} --new_mesh_dir {new_mesh_dir} --new_model_dir {new_model_dir} --model_tags {model_tags} --output_dir {output_dir}"

s = Slurm("interp", {"partition": "skx-normal",
                     "nodes": 10, "ntasks": 324, "time": "00:60:00", "account": "TG-EAR140030"})

s.run(f"date; {command1}; date;")
