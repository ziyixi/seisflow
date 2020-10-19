.. _handle_bin_files:

Handle bin files
==================================================================

The bin files are a type of file format in Specfem3D-globe to store the SEM (Spectral element method) mesh information. It can be 
divided into two parts:

*   **the "solver_data" files:** These type of files are named as "proc000000_reg1_solver_data.bin" in the "DATABASES_MPI" directory
    in the Specfem3D-globe simulation directory. It stores the location information of the SEM elements and the GLL points. (
    Gauss-Lobatto Legendre points) 
*   **the physical parameter files:** These files are usually named as "proc000440_reg1_beta_kernel.bin" for the "beta_kernel". The "solver_data" 
    files map the location of the elements and the GLL points into a 1D array, and the physical parameter is also a 1D array with the same size.
    
The operations performed on the bin files usually include:

*   Perform the mathematical operations on the bin files. (eg: sum up two kernel files)
*   Interpolate the bin files to the regular grid.
*   Interpolate the physical parameters of one mesh to another mesh.
*   Generate the mesh files for the Point Spread Function test (Describe in another section.)