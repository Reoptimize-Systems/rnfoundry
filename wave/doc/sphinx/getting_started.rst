Getting Started
***************

As sated previously, |TNshort| is is based in Matlab code [#f1]_, 
which can also run in Octave. |TNshort| requires the free MBDyn 
multibody modelling package to describe and simulate the motion of 
the system. To understand what background knowledge is required to 
use the toolbox see the Section :ref:`required-knowledge`.

.. [#f1] The |TNshort| functions are in pure Matlab, but the MBDyn package uses a function written in C++ to communicate with the MBDyn multibody modelling program.


The general workflow of |TNShort| is to first obtain the geometry of 
the device to be simulated and to use this to generate hydrodynamic 
data using a BEM solver. The toolbox provides an interface to the 
Nemoh BEM solver to assist with this (see :ref:`nemoh-interface` for 
a guide to generating this data using Nemoh). Generating the data 
using any other software, such as WAMIT is beyond the scope of this 
document. Once data has been generated in the normal format of the 
BEM solver it can usually be converted to a format suitable for use 
in |TNShort| using the `BEMIO` function.

Code Organisation
=================

The original WEC-Sim requires that all code and the Simulink models 
of the system be placed in one directory and that the `wecSim` 
function be run from that directory. This is not the case for 
|TNShort|, as the entire problem is defined using Matlab code. The 
only non code files required are the geometry and BEM output files 
(and the BEM data can be saved to a normal .mat file). However, the 
authors of |TNShort| do recommend a particular organisation which 
can be helpful and which uses the Matlab concept of 
:ref:`matlab-packages`. If you are not familiar with `packages` as 
a  way of organising code, you should first read this section, and 
the Matlab documentation
