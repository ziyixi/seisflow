.. _source_inversion:

Source inversion on stampede2
==================================================================

In Seisflow, we perform the source inversion using **the Adjoint Method** based on Specfem3D-globe. When we do the seismic tomography using FWI, 
it's usually important to invert both the source and structure, as they will all influence the waveform.

We perform the source inversion for the location and the focal mechanism, as well as the event time. Following the method in the paper `Adjoint centroid-moment tensor inversions <https://academic.oup.com/gji/article/186/1/264/698300>`__.
We divided the workflow into the following parts.

1.  Do the forward simulation of the current cmt solution and calculate the synthetics and the wave field.
#.  Calculate the misfit and the adjoint sources.
#.  Based on the method in adjoint centroid-moment tensor inversions, do the adjoint simulation and calculate the gradient for 
    each component of the focal mechanism and the location. (Note: The use of the development branch of Specfem3D-globe is needed as 
    in older versions, there is a bug that influences the gradient.)
#.  We take the gradients into two groups, one group contains the focal mechanism and another group contains the locations. In each group, 
    we use one common step length as the gradients are comparable within the group. For the two step lengths of the groups, we use a scaling 
    factor proposed in the method paper. 
#.  After getting the gradient, we use the steepest descent method to generate a perturbed cmt solution file of a fixed step length, and use that to calculate the 
    synthetics of the perturbed source.
#.  Calculate the synthetics of the perturbed source, if we assume the synthetics is linear to the step length we calculate the perturbed source, we can approximate the 
    synthetics of different step lengths and do the line search. Apart from the step length of the focal mechanism and the locations, we will also invert for another parameter
    named event time shift. Since we can just shift the synthetics to simulate the shift of the event time, we can just do the grid search to find the optimal step length 
    and the event time shift. This process is kind of time-consuming and we use the 
    method of the `Bayesian Optimization <https://github.com/fmfn/BayesianOptimization>`__ to do the grid search.

Prepare for the inversion
---------------------------------------

We also need the windows, the data-info and the processed data which are the same as the structure inversion. You can refer to the structure inversion to see how to prepare these files.
In the practice, these files are the same so we might just copied them to the root directory of the source inversion.

Execute the inversion script
----------------------------------

Here we give an example of the source inversion shell script which calls ``seisflow.scripts.xsede.xsede_perform_source_inversion``.

.. literalinclude :: ../_static/code/source_inversion.sh
   :language: bash

So everything is very similar to the structure inversion, but there are some flags that need to be explained:

*   **alpha_range:** To decide an appropriate search range of the step length, it's better to use the default range in the script. 
    As the meaning of the step length in the source inversion is difficult to explain (related to the magnitude of the gradient), it's better to 
    have a test when it's applied to the real inversion.
*   **tau_range:** Apart from the focal mechanism, the locations, and the event time, we can also invert for the half duration. However, according to my
    test, the misfit is almost insensitive to the half duration. According to the paper of the global cmt solution, the source half duration is simply determined
    by the magnitude of the source, so we should avoid inverting it.