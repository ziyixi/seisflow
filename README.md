# Seisflow

## Introduction

Seisflow is a package to perform the full seismic waveform inversion. It's designed to have the following features:
+ Integrate with the [Specfem3d globe](https://github.com/geodynamics/specfem3d_globe) package.
+ Be able to reuse the code and seprate the modules that are independent. (todo: use dynamic module loading for customization).
+ Use [ASDF](https://seismic-data.org/) to store and process the seismic data.
+ Integrate with Slurm, be able to submit the job automatically and correctly handle the job dependency. 

The package is provided with no warranty, I will do my best to keep the code usable and correct, but the user should also check the code and be responsible to their own results.

## Install
The code is avaliable from [Github](https://github.com/ziyixi/seisflow). Since this package is still under development, It's better to clone the repository directlly and put it into your working directory, 
i.e.:
```bash
cd [Your working directory]
git clone https://github.com/ziyixi/seisflow repo
ln -s repo/seisflow seisflow
```
So in the working directory, we can directly import the seisflow package and run the scripts using `python -m`.

If the package has been downloaded somewhere else, it's also ok to add the seisflow path to the python path:
```python
import sys
sys.path.append("[seisflow directory path]")
```

## Code structure
All the functions should be accessed from calling a script in the scripts directory, and these scripts could be called from the working directory, for example:
```bash
python -m seisflow.scripts.build_structure --help
```
And `seisflow.scripts.build_structure` is corresponding to `seisflow/scripts/build_structure.py`.

Below is the brief description for each subdirectory:
+ **asdf**: functions to handle asdf related problems, such as processing the data and convert sac files to the asdf format.
+ **data**: this directory is designed to store some data, for example, cmpaz_segment.txt is the file storing the orientation information of seismic stations from Chinese Array. 
+ **download**: scripts to download the data.
+ **julia**: directory to store the Julia scripts. To run the scripts, firstly you need to compile the fortran dependency:
    ```bash
cd seisflow/seisflow/julia/specfem_gll.jl/src/fortran/
julia compile.jl
    ```
    And later you can run scripts in `seisflow/seisflow/julia/scripts/`
+ **plot**: The scripts related to plotting maps or the figures.
+ **scripts**: The mostly used scripts.
+ **slurm**: A wrap of the package [slurmpy](https://github.com/brentp/slurmpy), also with bug fixed.
+ **tasks**: The functions that are called from scripts.
+ **utils**: Some utils functions.

## Workflow
The package is designed to automize as much as possibe but keep the flexbility. Below is a tutorial on using this package to perform the full seismic waveform inversion.

### Step1: prepare the data.
There are some scripts to convert seed files and sac files to the asdf format:
+ **seisflow/scripts/sac2asdf_hinet.py**
+ **seisflow/scripts/sac2asdf.py**
+ **seisflow/scripts/seed2asdf.py**
+ **seisflow/scripts/mpi_sac2asdf_hinet.py**
+ **seisflow/scripts/mpi_sac2asdf.py**

Generally speaking, the scripts with the name "mpi_" at the beginning is a job that should be run in parallel. (with `mpirun -np`). This scripts will never handle all the cases and
you may need to modify them to meet your requirement. A good reference will be the [pyasdf document](http://seismicdata.github.io/pyasdf/).

### Step2: Process the data.
From step1, you may get different types of the asdf files. The difference is how the response files are provided. If they are in the PZ format, you may have to refer to mpi\_sac2asdf\_hinet.py, which will simply convert the sac files to the asdf format without the repsonse information. If they are 
provided in RESP format, which means it can be stored in Obspy and further pyasdf, you should refer to mpi_sac2asdf.py, which will have the Stationxml information stored. Thus we don't have to provide the 
response information during processing the data for this kind of asdf file.

There are several processing scripts:
+ **seisflow/scripts/process_data_pz.py**: Process the asdf file with PZ.
+ **seisflow/scripts/process_data.py**: Process the asdf file that has resp information in it.
+ **seisflow/scripts/process_sync.py**: Process the sync asdf files. (just do the frequency taper with no convolution)
+ **seisflow/scripts/xsede_process_sync.py**: Process the sync asdf files in stampede2.
+ **seisflow/scripts/icer_process_sync.py**: Process the sync asdf files in ICER.
+ **seisflow/scripts/mpi_simplify_data.py**: Remove the station in an asdf file that has no data in it.

### Step3: Extract data information and build up windows.
Since there are some data information that are widely used, such as the theoritical first phase arrival time, we have to extract the data information from the asdf files and store them in some place.
Calling the script `seisflow/scripts/mpi_extract_data_info.py` will do this trick. Also you can modify the function `extract_data_info`. The output python puckle files are stored as one pkl file for one category, which means you can add new categories easilly.

The package has provided a script `seisflow/scripts/mpi_tao_2018_ggg_windows.py` to build windows with fixed length. All the windows should be inherited from the `Window` class in `seisflow/tasks/windows/window.py`.
It's easy to write a script to convert the outputed windows from Flexwin to the windows following the same structure as described in `seisflow/scripts/mpi_tao_2018_ggg_windows.py`.

### Step4: Prepare the synthetics.
I will write a blog later to talk about how to enable ASDF in Specfem3D globe to output the asdf syncthetics. Assume we have prepared and compiled the reference specfem directory, whose bin subdirectory will be used for each single event. You can use `seisflow/scripts/build_structure.py` to build up the inversion structure and `seisflow/scripts/xsede_run_multiple_forward_jobs.py` to submit the simulation slurm script in stampede2.

After finishing calculation, you can use `seisflow/scripts/collect_src_frechet.py` to collect all the synthetics to an output directory.

### Step5: Perform the source inversion.
To perform the source inversion following [Tao's method](seisflow/scripts/collect_src_frechet.py), user can refer to `seisflow/scripts/xsede_perform_source_inversion.py`. One important thing is that **ALL THE PATHS IN THE SCRIPT SHOULD BE THE ABSOLUTE PATH**. (to avoid relative path related bugs.) All the flags in the script is self-explainable, and there are some flags that may need some explanations:
+ **inversion_directory**: all the inversion will be performed within this directory.
+ **raw\_sync\_directory**: the calculation of the adjoint source should have the same time steps as being set in the synthetics directly outputed from Specfem. Thus this directory will contain all the "raw" synthetics without processing.

The recommended way of performing the source inversion is:
1. make a directory as the base directory.
2. clone the package following the installation part.
3. prepare all the required directories in `xsede_perform_source_inversion.py`. It's better to either make a soft link or copy the directories to the base directory.
4. run the script to submit the job.

One thing to note is that the job is for one single iteration. Usually after two or three iterations, the CMT solution will converge. Since we have to run at least two iterations, the recommended way will be:
1. Submit the first job.
2. At the same time, submit the second job.
3. Cancel the second job.
4. In the slurm-scripts directory, there will be the job script for the second iteration. Thus based on Slurm's dependency, we can submit the job using `sbatch --dependency:afterok:[jobid1] slurm-scripts/[job script 2]`. Note the slurm script must be submitted in the base directory but not in the slurm-script directory as the slurm script is using relative path.

### Step6: Perform the structure inversion.
User can refer to `seisflow/scripts/structure_inversion/xsede_perform_structure_inversion.py` to perform the structure inversion, all the flags should also be self-explainable. The setting is very similar to the source inversion part. And also there are something that should be kept in mind:
1. The script to calculate the adjoint sources are different. Since performing structure inversion will weight for all the events.
2. The script will need extra package [sem_utils](https://github.com/taotaokai/sem_utils) to do kernel smoothing.

## Summary
In the future, I will write some blogs to discuss the architectures of designing the source and structure inversion directories, which is very similar to `React.js`'s JSX way of writing code. For the future development, maybe it's better to seprate the logic apart from the implementation. 