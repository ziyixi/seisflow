.. _specfem:

Configuring Specfem3D-globe
=================================

Why ASDF?
---------------

To use the `ASDF <https://seismic-data.org//>`__ feature of Specfem3D-globe, an additional compilation configuration is 
required. As for the reasons why we prefer ASDF format but not the traditional SAC format, we quote its benefits from the ASDF website:

1. The amount of seismic data available for analysis worldwide is rapidly growing. Seismic arrays, such as USArray and ChinaArray, give access to datasets on the terabyte scale that are not suited for existing seismic data formats.
#.  Disk space is rapidly growing and data organization should improve such that the different types of seismic data (waveforms, receivers, earthquakes, adjoint sources, cross correlations, etc.) can be easily exchanged among the community under one container.
#.  Modern workflows in seismology use supercomputers and the number of files is an I/O bottleneck. The performance of these workflows would be increased if the data was stored by combining all time series into one file and taking advantage of parallel processing capabilities.
#.  New methods, such as ambient-noise seismology, should not be limited by data formats that were developed for other applications in seismology. Also, seismologists often ignore standards because adherence increases development time. An adaptable seismic data format with an open, modular design will be able to evolve and handle future advances in seismology.
#.  Reproducibility is a goal in science and seismology has yet to develop a standardized way of storing provenance in the current seismic data formats. We introduce a format that contains flexible provenance that lets the user know where the data comes from and what has been done to it.

In my experience, using ASDF will generate less files in the simulation, and since in HPC it usually limits the quota of file numbers, 
using one ASDF file instead of lots of SAC files will help with this issue. Additionally, it's easy to back up the files and transfer between
different systems. By using the Pyasdf package (Thanks to Lion Krischer), we can also process the ASDF files in parallel, which is a more elegant 
and easier way to code.

Enable the ASDF feature of Specfem3D-globe
------------------------------------------------------------

Here we take the development version of Specfem3D-globe as an example. Other versions later than v7.0.0 will be slightly different.

Compile ASDF
>>>>>>>>>>>>>>>>>>>>>

The Specfem3D-globe package is linked with a static library named libasdf.a . To get the static library, we have to compile from the source code.
The detail could be referred from `the ASDF document <https://github.com/SeismicData/asdf-library>`__ 

One thing to note is that here we should use the parallel version of HDF5. On some HPC systems such as stampede2, the module system has provided 
the parallel HDF5, we can directly load it. However, in some other systems, we might have to compile the parallel version by ourself. The detail of
compiling parallel HDF5 could refer to `The document of HDF5 <https://support.hdfgroup.org/HDF5/PHDF5/>`__

Edit Parfile
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Before compilation, the first step is to edit the Parfile in the DATA directory. One feature of Specfem3D-globe is that the editing of 
the Parfile might result in a recompilation. It's because some contents in Parfile will be included in some generated ``.h`` files and these 
files will be included in the source code in the compilation. The configuration of Parfile could refer to `the Document <https://geodynamics.org/cig/software/specfem3d_globe/>`__.
We should let ``OUTPUT_SEISMOS_ASDF`` and ``READ_ADJSRC_ASDF`` to be true.

Compile Specfem3D-globe
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

The detail of compiling Specfem3D-globe could refer to the `Document <https://geodynamics.org/cig/software/specfem3d_globe/>`__. The main difference between 
is that we will need to enable ASDF in configuration. Here we take the example in stampede2::

    # config to use the asdf files
    HDF5=/opt/apps/intel18/impi18_0/phdf5/1.8.16/x86_64/lib
    ASDF=/work/05880/tg851791/stampede2/asdf-library-1.0.0/asdf/lib

    ./configure FC=ifort CC=icc CXX=icpc MPIFC=mpif90 ASDF_LIBS="-L${ASDF} -L${HDF5} -lpthread"  --with-asdf

Here I install the ASDF library in ``/work/05880/tg851791/stampede2/asdf-library-1.0.0/asdf/lib``, and I am loading phdf5-1.8.16. And 
The Intel compilers have also been used. The content in ``ASDF_LIBS`` is a trick as in the compilation, there might be the errors like some 
libraries are not found, in which cases we should use ``-L`` or ``-l`` to provide the  path of the libraries.
