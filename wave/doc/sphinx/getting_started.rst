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


Code and Data Organisation
==========================

The original WEC-Sim requires that all code and the Simulink models 
of the system be placed in one directory and that the `wecSim` 
function be run from that directory. This is not the case for 
|TNShort|, as the entire problem is defined using Matlab code. The 
only non code files required are the geometry and BEM output files 
(and the BEM data can be saved to a normal .mat file). However, the 
authors of |TNShort| do recommend a particular organisation which 
can be helpful.

Code Organisation
-----------------

The code organisation method described uses the Matlab concept of 
packages. If you are not familiar with packages as a  way of 
organising code, you should first read the section 
:ref:`required-knowledge-matlab-packages`.

However, in summary, packages are created by placing files in a 
directory whose name begins with a ``+`` character. For example, one 
might create a directory named ``+mypackage``. Any script functions 
put in this directory can be used from the Matlab prompt using the 
syntax::

   mypackage.function_name 
   
where ``function_name`` is the name of the script or function file 
in the package directory. This allows you to have the same function 
and script names within different packages without clashing (also 
called shadowing). Organising things this way makes it easier to 
make new designs/packages of scripts and functions from existing 
ones by just copying the package directory to a new directory (also 
starting with a ``+`` symbol). We will use this method of 
organisation in all subsequent examples.

The package, or function files used to create the model can reside 
anywhere provided they are known by Matlab (i.e. they are on the 
Matlab path, see :ref:`required-knowledge-matlab-packages` for more 
information). More information on what files are actually needed to 
generate the model will be shown below through an example.

Data Organisation
-----------------

In contrast to the code, a certain structure is required for the 
organisation of the input data files used by |TNshort|. Data is 
organised in a case directory which must have two subdirectories, 
"geometry" and "hydroData".

| project
| ├── geometry
| ├── hydroData


geometry
   Contains `STL`_ files for each of the (hydrodynamically 
   interacting) bodies in the system. STL is a mesh format. These 
   files are used for visualisation and some hydrodynamic 
   calculations during a simulation, e.g. nonlinear buoyancy forces.

.. _STL: https://en.wikipedia.org/wiki/STL_(file_format)

hydroData
   Contains hydrodynamic data files, which can be of various types, 
   and may include subdirectories. Processed hydrodynamic data files 
   are also stored here. Nemoh, WAMIT and AQUA BEM solver output 
   files can be processed into the required inputs for |TNShort| 
   using the BEMIO functions. At a minimum, |TNShort| requires that 
   there is either one HDF5 (.h5) formt file (which can be produced 
   using the Write_H5 function), or a set of .mat files, one for 
   each body (which can be produced using the 
   ``wsim.bemio.write_hydrobody_mat_files`` function). Currently 
   only the set of mat files format is possible when using Octave 
   due to missing functionality to load the HDF5 format files 
   (specifically the h5read function). 
   
As explained above, it is not necessary for the code files to be 
located in the same place as the data files. However, similarly, 
there is no reason they may *not* be located in the project folder, 
and furthermore, the case directory can also be a package folder such 
that you end up with a directory tree something like:

| +project_name
| ├── geometry
| ├── hydroData
|    ├── body_id_0.mat
|    ├── body_id_1.mat
|    ├── body_id_3.mat
| ├── generate_hydrodata.m
| ├── run.m
| ├── make_multibody_system.m

Where ``project_name`` is the Matlab packge name, so the functions 
within it are called like ``project_name.generate_hydrodata``, 
``project_name.run`` and ``project_name.make_multibody_system`` 
within Matlab.


Example 1: The RM3 Two Body Point Absorber
==========================================

This section describes the application of the WEC-Sim code to model 
the Reference Model 3 (RM3) two-body point absorber WEC. This 
example application can be found in the |TNshort| examples 
directory. In this example, we will start from having only an 
existing geometry file, to simulating the entire system.

Device Geometry
---------------

The RM3 two-body point absorber WEC has been characterized both 
numerically and experimentally as a result of the DOE-funded 
Reference Model Project. The details and outcomes of this study can 
be found `here`__. The RM3 is a two-body point absorber consisting 
of a float and a reaction plate. Full-scale dimensions of the RM3 
and its mass properties are shown below.

.. __: http://energy.sandia.gov/energy/renewable-energy/water-power/technology-development/reference-model-project-rmp/

.. image:: /images/RM3_Geom.png

|

+--------------------------------------------------------------+
| Float Full Scale Properties                                  |
+--------+--------------+--------------------------------------+
| CG (m) | Mass (tonne) | Moment of Inertia (kg-m^2)           |
+========+==============+============+============+============+
| 0.0    |              | 20'907'301 |            |            |
+--------+              +------------+------------+------------+
| 0.0    |  727.0       |            | 21'306'091 | 4305       |
+--------+              +------------+------------+------------+
| -0.72  |              |            | 4305       | 37'085'481 |
+--------+--------------+------------+------------+------------+

|

+--------------------------------------------------------------+
| Plate Full Scale Properties                                  |
+--------+--------------+--------------------------------------+
| CG (m) | Mass (tonne) | Moment of Inertia (kg-m^2)           |
+========+==============+============+============+============+
| 0.0    |              | 94'419'615 |            |            |
+--------+              +------------+------------+------------+
| 0.0    |  878.30      |            | 94'407'091 | 217'593    |
+--------+              +------------+------------+------------+
| -21.29 |              |            | 217'593    | 28'542'225 |
+--------+--------------+------------+------------+------------+


Generating Hydrodynamic Data
============================

The first step in modelling the system is to generate the 
hydrodynamic data files using a BEM solver such as Nemoh, WAMIT or 
AQUA. This will be demonstrated in this example using Nemoh. At this 
point it may be worth reading the section :ref:`nemoh-interface` 
which intrduces the Matlab based preprocessor which has been 
developed to help with this with simple examples.

To use Nemoh, you must first generate surface meshes of the bodies 
you wish to simulate. 
