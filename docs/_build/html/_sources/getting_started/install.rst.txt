.. _install:

Installing
==========

The environments
--------------------

You'll need **Python 3.6 or greater** to use seisflow (as f-string is widely used in the package).

Using `Anaconda <https://www.anaconda.com/distribution>`__ is recommended as it could ensure all the dependencies to be installed properly. And it's also the 
development environment I'm currently using.

The **parallel version** of `Pyasdf <https://seismicdata.github.io/pyasdf/>`__ is required when the data or synthetics processing is included (such as doing the structure inversion).
It could not be directly downloaded through anaconda, and the detail could refer to the documentation of Pyasdf.

The `Specfem3D-globe <https://geodynamics.org/cig/software/specfem3d_globe/>`__ should **support ASDF** in compilation. The details will be discussed in the next section.

Requirements
---------------------

.. include:: ../../requirements.txt
   :literal:

Install seisflow
----------------------

Since the package is under development and aimed for the code modification by users, the recommended usage is to directly clone this package into a directory and use 
it there::

    cd <your directory>

And then we download the package, but rename it to `repo` since the seisflow directory inside `repo` is the python package::

    git clone https://github.com/ziyixi?tab=repositories repo 

Here we could make a soft link::

    ln -s repo/seisflow .

Test Installation
--------------------------

To test the installation, in the directory you clone the package, you can run::

    python -m seisflow.check_installing

It checks whether you have installed all the required packages. If everything is OK, it will output::

    All your dependencies have been installed.