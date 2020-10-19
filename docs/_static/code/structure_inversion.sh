#!/bin/bash

PY=/work/05880/tg851791/stampede2/anaconda3/envs/seismology/bin/python # the python path
BASE=/scratch/05880/tg851791/big_region_structure_inversion # the base directory
ITER=25 # iteration number
CURRENT=${BASE}/iter${ITER} # the root directory

base_directory=${CURRENT}/inversion # all the simulation will be performed inside it
cmts_directory=${CURRENT}/cmts # the cmt solution directory
ref_directory=${CURRENT}/ref # the reference specfem3D-globe directory
windows_directory=${CURRENT}/windows # the windows directory
data_asdf_directory=${CURRENT}/data_asdf # the processed data directory
data_info_directory=${CURRENT}/data_info # the data info directory
last_step_model_update_directory=${CURRENT}/last_step_model_update # the directory storing the gradient for the last iteration
stations_path=${CURRENT}/stations # the stations path
sem_utils_directory=${CURRENT}/sem_utils # the directory of sem_utils
source_mask_directory=${CURRENT}/source_mask # the directory of the external source mask

n_total=142 # total number of events
n_each=13 # for each iteration of the modeling, how many events we are running
n_iter=11 # the total number of iterations we are running
nproc=441 # for each event, how many MPI processes we are using
n_node=120 # the total number of nodes
partition=skx-normal # the partition in stampede2
time_forward=20:00:00 # the time of the forward modeling
account=TG-EAR140030 # the account of the project
n_node_process_kernel=10 # the number of nodes to process the kernel
time_process_kernel=02:00:00 # the time of processing the kernel
time_run_perturbation=08:00:00 # the time of calculating the perturbed model and doing the line search
periods=8,40/20,120 # the filtering time periods, / seprate the body and surface waves
waveform_length=1800 # the waveform length in seconds
sampling_rate=10 # the sampling rate
taper_tmin_tmaxs=1,50/1,150 # the tapering time, / seprate the body and surface waves
sigma_h=25 # the horizontal raidus of gaussian smoothing
sigma_v=10 # the vertical raidus of gaussian smoothing
search_range=0,0.05 # the line search range, from 0 to 5%
search_step=0.002 # the line search steo length

# * run the inversion script
${PY} -m seisflow.scripts.xsede.xsede_perform_structure_inversion --base_directory ${base_directory} --cmts_directory ${cmts_directory} --ref_directory ${ref_directory} --windows_directory ${windows_directory} --data_asdf_directory ${data_asdf_directory} --data_info_directory ${data_info_directory} --last_step_model_update_directory ${last_step_model_update_directory} --stations_path ${stations_path} --sem_utils_directory ${sem_utils_directory} --source_mask_directory ${source_mask_directory} --n_total ${n_total} --n_each ${n_each} --n_iter ${n_iter} --nproc ${nproc} --n_node ${n_node} --partition ${partition} --time_forward ${time_forward} --account ${account} --n_node_process_kernel ${n_node_process_kernel} --time_process_kernel ${time_process_kernel} --time_run_perturbation ${time_run_perturbation} --periods ${periods} --waveform_length ${waveform_length} --sampling_rate ${sampling_rate} --taper_tmin_tmaxs ${taper_tmin_tmaxs} --sigma_h ${sigma_h} --sigma_v ${sigma_v} --search_range ${search_range} --search_step ${search_step}