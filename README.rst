Seisflow
==========

Seisflow is an FWI (Full Seismic Waveform Inversion) workflow package based on `Specfem3D-globe <https://geodynamics.org/cig/software/specfem3d_globe/>`__
and `Pyasdf <https://seismicdata.github.io/pyasdf/>`__.

`Documentation (development version) <https://www.pygmt.org/dev>`__ |
`Contact <https://forum.generic-mapping-tools.org>`__

.. placeholder-for-doc-index

About
-------------

🚨 **This package is still undergoing rapid development.** 🚨

Seisflow is a package for building up the workflow in FWI and has provided lots of useful scripts that are commonly used in the inversion which are highly embedded into the 
full workflow. 

Seisflow has mainly provided functions as:

* Handle the mesh files of Specfem3D-globe.
* A workflow to perform the source inversion based on the adjoint method.
* A workflow to perform the structure inversion based on the adjoint method.
* Scripts to handle ASDF files. (a a conversion with SAC, processing, etc.)
* Scripts to download the data. (from hinet, etc.)
* Scripts to do the Point Spread Function (PSF) test to get the model resolution.
* Plotting Scripts for the generated model based on `pygmt <https://www.pygmt.org/>`__

Designs
----------

Due to the complexity of the FWI, it's hard to design a package that meets the needs of everyone. When I develop the package, I intend to split the whole package into small parts.
And in each part, it might be a python file that you can modify the content. By exposing the same API, there is no need to concern how other parts are designed. For example, if you 
want to use a new misfit function, you can just modify one function to calculate the misfit windows, and it can be assured that this function is the only part to calculate the misfit and other
functions will need to import it.

License
-------

Seisflow is free software: you can redistribute it and/or modify it under the terms of
the **BSD 3-clause License**. A copy of this license is provided in
`LICENSE <https://github.com/ziyixi/seisflow/blob/master/LICENSE>`__.

Support
-------

The development of Seisflow is funded by
`NSF grant EAR-1802247 <https://www.nsf.gov/awardsearch/showAward?AWD_ID=1802247>`__.