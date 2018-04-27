Getting Started
***************

As sated previously, |TNshort| is is based in Matlab code [#f1]_, 
which can also run in Octave. |TNshort| requires the free MBDyn 
multibody modelling package to describe and simulate the motion of 
the system. To understand what background knowledge is required to 
use the toolbox see the Section :ref:`required-knowledge`.



The general workflow of |TNShort| is to first obtain the geometry of 
the device to be simulated and to use this to generate hydrodynamic 
data using a BEM solver. The toolbox provides an interface to the 
Nemoh BEM solver to assist with this (see :ref:`nemoh-interface` for 
a guide to generating this data using Nemoh). Generating the data 
using any other software, such as WAMIT is beyond the scope of this 
document. Once data has been generated in the normal format of the 
BEM solver it can usually be converted to a format suitable for use 
in |TNShort| using the `BEMIO` function.


.. [#f1] The |TNshort| functions are in pure Matlab, but the MBDyn 
   package uses a mex function written in C++ to communicate with 
   the MBDyn multibody modelling program. Uers don't need to 
   directly interact with C++ code to use the toolbox.


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
:ref:`required-knowledge-matlab-packages`. If you are not familiar 
with packages as a  way of organising code, you should first read 
this section, and the `Matlab documentation`__.

.. __: https://uk.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html

Packages are created by placing files in a directory whose name 
begins with a ``+`` character. For example, one might create a 
directory named ``+mypackage``. Any script functions put in this
directory can be used from the Matlab prompt using the syntax::

   mypackage.function_name 
   
where ``function_name`` is the name of the fscript or function file 
in the package directory. This allows you to have the same function 
and script names within different packages. Organising things this 
way makes it easier to make new designs/packages of scripts and 
functions from existing ones by just copying the package directory 
to a new directory (also starting with a ``+`` symbol). We will use 
this method of organisation in all subsequent examples. 

Example 1: The RM3 Point Absorber
=================================



