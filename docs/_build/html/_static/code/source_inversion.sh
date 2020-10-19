#!/bin/bash

PY=/work/05880/tg851791/stampede2/anaconda3/envs/seismology/bin/python # the python path
BASE=/scratch/05880/tg851791/big_region_structure_inversion # the base directory
CURRENT=${BASE}/source_inversion_after_iter20 # the root directory

inversion_directory=${CURRENT}/inversion # the inversion directory
cmtfiles_directory=${CURRENT}/cmts # the cmtfiles directory
ref_directory=${CURRENT}/ref # the specfem reference directory
data_directory=${CURRENT}/data_asdf # the processed data directory
windows_directory=${CURRENT}/windows # the windows directory
data_info_directory=${CURRENT}/data_info # the datainfo directory
stations_path=${CURRENT}/stations # the stations path
raw_sync_directory=${CURRENT}/raw_sync # the raw sync directory, we can just generate a set of ".h5" files for all the events as if they are raw synthetic files.

n_total=142 # total number of events
n_each=13 # number of events for each iteration in specfem simulation
n_iter=11 # number of iterations in specfem simulation
nproc=441 # the number of processes for each specfem simulation
n_node=120 # the total number of nodes to use
ntasks=5733 # the total number of mpi processes to use
partition=skx-normal # the partion to use
simulation_time_step1=24:00:00 # the simulation time for step 1 (forward+adjoint+forward for perturbed model)
account=TG-EAR140030 # the account to use
n_node_line_search=10 # the total number of nodes to use in step 2 (grid search for the step length and the event time shift)
ntasks_line_search=142 # the total number of processes to use in step 2
partition_line_search=skx-normal # the partion to use in step 2
simulation_time_step2=02:00:00 # the simulation time for step 2
waveform_length=1800 # the waveform length
taper_tmin_tmaxs=1,50/1,150 # the taper time bands: minp1,maxp1/minp2,maxp2/...
periods=8,40/20,120 # periods in filtering: minp1,maxp1/minp2,maxp2/...
sampling_rate=10 # the sampling rate to use
alpha_range=-10,10 # the line search range for alpha (step length for the focal mechanism and the location)
t0_range=-3,3 # the line search range for t0
tau_range=fixed # the line search range for tau, use 'fixed' if no change

${PY} -m seisflow.scripts.xsede.xsede_perform_source_inversion --iter_number ${iter_number} --py ${PY} --n_total ${n_total} --n_each ${n_each} --n_iter ${n_iter} --nproc ${nproc} --n_node ${n_node} --ntasks ${ntasks} --partition ${partition} --simulation_time_step1 ${simulation_time_step1} --account ${account} --n_node_line_search ${n_node_line_search} --ntasks_line_search ${ntasks_line_search} --partition_line_search ${partition_line_search} --simulation_time_step2 ${simulation_time_step2} --inversion_directory ${inversion_directory} --cmtfiles_directory ${cmtfiles_directory} --ref_directory ${ref_directory} --data_directory ${data_directory} --windows_directory ${windows_directory} --data_info_directory ${data_info_directory} --stations_path ${stations_path} --raw_sync_directory ${raw_sync_directory} --waveform_length ${waveform_length} --taper_tmin_tmaxs ${taper_tmin_tmaxs} --periods ${periods} --sampling_rate ${sampling_rate} --alpha_range ${alpha_range} --t0_range ${t0_range} --tau_range ${tau_range}
