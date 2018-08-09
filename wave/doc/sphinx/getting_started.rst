Getting Started
***************

As stated previously, |TNshort| is is based in Matlab code [#f1]_,
which can also run in Octave. |TNshort| also requires the free MBDyn
multibody modelling package to describe and simulate the motion of
the system (it is likely you will have received this alongside
|TNShort|). The MBDyn Matlab Toolbox also has its own manual with
simple examples demonstrating the toolbox's capabilities. As
|TNShort| makes heavy use of this toolbox, it is advisable to first
examine some of the examples in this manual to get a better
understanding of how the pre-processor works before attempting to
develop a system for |TNShort| as an in-depth discussion of the
operation and organisation of the MBDyn preprocessor will not be
provided here.

To understand what other background knowledge is required to
use |TNshort| see the Section :ref:`required-knowledge`.

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
directory whose name begins with a `+` character. For example, one
might create a directory named `+mypackage`. Any script functions
put in this directory can be used from the Matlab prompt using the
syntax::

   mypackage.function_name

where ``function_name`` is the name of the script or function file
in the package directory. This allows you to have the same function
and script names within different packages without clashing (also
called shadowing). Organising things this way makes it easier to
make new designs/packages of scripts and functions from existing
ones by just copying the package directory to a new directory (also
starting with a `+` symbol). We will use this method of
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
   calculations during a simulation, e.g. non-linear buoyancy forces.

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

Where `project_name` is the Matlab packge name, so the functions
within it are called like

::

   project_name.generate_hydrodata

and

::

   project_name.run

and

::

   project_name.make_multibody_system

within Matlab.


Example: The RM3 Two Body Point Absorber
==========================================

This section describes the application of the WEC-Sim code to model
the Reference Model 3 (RM3) two-body point absorber WEC. This
example application can be found in the |TNshort| examples directory
in the +example_rm3 subdirectory. In this example, we will start
from having only an existing geometry file, to simulating the entire
system.

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
which introduces the Matlab based preprocessor which has been
developed to help with this with simple examples.

To use Nemoh, you must first generate surface meshes of the bodies
you wish to simulate. Ideally these must be quadrilateral meshes,
although triangular meshes can be treated as degenerate
quadrilaterals. These meshes must be cut off at the mean free
surface level of the water, so they will likely be different from
the STL meshes used for visualisation (unless all bodies are not
surface piercing). The mesh for the RM3 system for Nemoh is shown in
:numref:`nemoh_rm3_mesh` and the mesh files are provided in the
geometry subdirectory. The Nemoh mesh input is defined using two
files, a .nmi and a .cal file.

.. _nemoh_rm3_mesh:
.. figure:: /images/nemoh_rm3_mesh.png
   :align: center

   Plot of both Nemoh input meshes for the RM3 example.

In this example, generating the hydrodynamic data is performed by
the `example_rm3.generate_hydrodata` function. We will now work
our way through this function to show how the process works (but
remember you should first read the section
:ref:`nemoh-interface` to make understanding this easier).

The start of this function just sets up from option parsing and
directories etc. which aren't important for understanding it's
actual operation:

.. literalinclude:: /examples/+example_rm3/generate_hydrodata.m
   :end-before: %% create the float

See the help for the `parse_pv_pairs` function to understand this
option parsing.

The next step is to create the Nemoh body objects for the float and
spar, using as input the mesh files for each:

.. literalinclude:: /examples/+example_rm3/generate_hydrodata.m
   :start-at: %% create the float
   :end-before: %% Create the nemoh simulation

The Nemoh bodies are now ready to be inserted into a Nemoh system
which will process their meshes and run the Nemoh calculations.

.. literalinclude:: /examples/+example_rm3/generate_hydrodata.m
   :start-at: %% Create the nemoh simulation
   :end-before: %% draw the course body mesh (this will be refined later)

The mesh can then be plotted, and it is this plot which was shown in
:numref:`nemoh_rm3_mesh`.

.. literalinclude:: /examples/+example_rm3/generate_hydrodata.m
   :start-at: %% draw the course body mesh (this will be refined later)
   :end-before: % write out the course mesh files for all bodies

Now the meshes can be processed. This involves writing out the Nemoh
mesher input files and calling `processMeshes` which runs the
Nemoh mesher on the the course mesh input files to produced a
refined mesh, and also performs some basic hydrodynamic and
geometrical calculations, with the results being put in the
hydroData directory.

.. literalinclude:: /examples/+example_rm3/generate_hydrodata.m
   :start-at: % write out the course mesh files for all bodies
   :end-before: %% Create hydro structure

In this case we have chosen to obtain results for 260 wave
frequencies. This will take quite a long time to run on a typical
desktop PC or laptop, on the order of several hours, and produce
several hundred MB of data in the form of text files. Generally you
only want to produce these results once.

The final part of the function deals with converting the output from
Nemoh into a form suitable for use in |TNShort|.

.. literalinclude:: /examples/+example_rm3/generate_hydrodata.m
   :start-at: %% Create hydro structure

The first step is to load the Nemo data files and convert them to a
standard format. This is achieved with the
`wsim.bemio.processnemoh` function. This function is part of a
package called bemio, which is within the wsim package. Note that
this is separate from the standard `BEMIO`_ functions from the
original WEC-Sim project. Copies of these are provided with
|TNShort| for convenience, any new function made specifically for
|TNShort| can be found in this `wsim.bemio` package.

.. _BEMIO: http://wec-sim.github.io/WEC-Sim/advanced_features.html#bemio

The output of `wsim.bemio.processnemoh` is a Matlab structure
containing hydrodynamic data for all the bodies in the Nemoh system.
This hydro structure can then undergo further processing as shown in
the function, depending on what types of simulation are desired to
be run in |TNShort|. Here we generate the wave radiation impulse
response functions for both the convolution integral calculation
method and the state-space form, and also the wave excitation
impulse response function.

Once all processing of the hydrodynamic data is complete, the data
must be converted to a form with can be used within |TNShort|, i.e.,
which can be used by the `wsim.hydroBody` class. This class will
be discussed in more detail in the following section, but for now it
is sufficient to know that the data must either be converted to on
HDF5 format file, or a set of .mat files, one for each
hydrodynamically interacting body. The first format is the same as
used by the original WEC-Sim project. The ability to use normal .mat
files directly has been added to |TNShort|, therefore only the .h5
file can be used with both tools. Compared to the output of Nemoh,
the .h5 and .mat files are quite small. The .h5 file will be on the
order of tens of MB (7.4MB in this case), while the .mat files will
be smaller still (two 1MB files in this case).

System Simulation
=================

Having generated the required hydrodynamic data to calculate the
wave interaction forces we are now ready to use this data in a time
domain simulation of the RM3 device. The time domain simulation and
system is defined in this case in one function,
`example_rm3.make_multibody_system` and a script
`example_rm3.run_wecsim`. Note that there is nothing special about
the names of these functions, you are free to organise the code in
any way you like, and name the scripts and functions any way you like.

Initial Simulation Setup
------------------------

We will start by examining `example_rm3.run_wecsim`, which in this
case is the main script which sets up and runs the simulation. The
first section of this script sets the general simulation and wave
parameters.

.. literalinclude:: /examples/+example_rm3/run_wecsim.m
   :end-before: %% Hydrodynamic body system

The simulation and wave settings are each stored using classes,
`wsim.simSettings` and `wsim.waveSettings`, where the settings
are the class properties. Most of the settings in this example are
self-explanatory, but more settings are available which may not be
as obvious. A full list of the possible simulation settings and
their descriptions are shown below:

+------------------+---------------------------------------------------------------------------------------------------------------------------+
| Property/Setting | Description                                                                                                               |
+==================+===========================================================================================================================+
| startTime        | Simulation start time (default = 0 s)                                                                                     |
+------------------+---------------------------------------------------------------------------------------------------------------------------+
| endTime          | Simulation end time (default = 500 s)                                                                                     |
+------------------+---------------------------------------------------------------------------------------------------------------------------+
| dt               | Simulation time step (default = 0.1 s)                                                                                    |
+------------------+---------------------------------------------------------------------------------------------------------------------------+
| dtFeNonlin       | Sample time at which to calculate nonlinear forces (default = dt)                                                         |
+------------------+---------------------------------------------------------------------------------------------------------------------------+
| dtCITime         | Sample time at which to calculate Convolution Integral (default = dt)                                                     |
+------------------+---------------------------------------------------------------------------------------------------------------------------+
| rampT            | Ramp time for wave forcing (default = 100 s)                                                                              |
+------------------+---------------------------------------------------------------------------------------------------------------------------+
| domainSize       | Size of free surface and seabed. This variable is only used for visualization (default = 200 m)                           |
+------------------+---------------------------------------------------------------------------------------------------------------------------+
| CITime           | Convolution integral time span (default = 60 s)                                                                           |
+------------------+---------------------------------------------------------------------------------------------------------------------------+
| ssCalc           | Option for convolution integral or state-space calculation: convolution integral->0, state-space->1, (default = 0)        |
+------------------+---------------------------------------------------------------------------------------------------------------------------+
| nlHydro          | Option for nonlinear hydrodynamics calculation: linear->'0', nonlinear->'1', (default = 0)                                |
+------------------+---------------------------------------------------------------------------------------------------------------------------+
| b2b              | Option for body-to-body interactions: off->false, on->true, (default = false)                                             |
+------------------+---------------------------------------------------------------------------------------------------------------------------+
| paraview         | Option for writing vtp files for paraview visualization.                                                                  |
+------------------+---------------------------------------------------------------------------------------------------------------------------+
| adjMassWeightFun | Weighting function for adjusting added mass term in the translational direction (default = 2)                             |
+------------------+---------------------------------------------------------------------------------------------------------------------------+
| mcrCaseFile      | mat file that contain a list of the multiple conditions runs with given conditions                                        |
+------------------+---------------------------------------------------------------------------------------------------------------------------+
| morrisonElement  | Option for Morrison Element calculation: Off->'0', On->'1', (default = 0)                                                 |
+------------------+---------------------------------------------------------------------------------------------------------------------------+
| rho              | Density of water (default = 1000 kg/m^3)                                                                                  |
+------------------+---------------------------------------------------------------------------------------------------------------------------+
| g                | Acceleration due to gravity (default = 9.81 m/s)                                                                          |
+------------------+---------------------------------------------------------------------------------------------------------------------------+

For the wave settings, examples of the different types of settings
are shown in the comments in this section of the script.


Creating the Hydrodynamic System
--------------------------------

The next step is to define the system of hydrodynamically
interacting bodies. This is done using the ``wsim.hydroBody`` class
and the ``wsim.hydroSystem`` class. The ``wsim.hydroSystem`` class
is essentially a container for a collection of hydroBody objects,
with one hydroBody object for each hydrodynamic body in the system.
We specifically refer to hydrodynamically interacting bodies a other
bodies in the system are defined elsewhere, as part of the multibody
system dynamics system. This will be explained in a later section.

The `wsim.hydroBody` class does all processing of the hydrodynamic
forces on a body. The `wsim.hydroSystem` organises these bodies,
inserts the correct position, velocity and acceleration information
into each body, and gets all of the forces applied to the bodies by
the wave action as the simulation proceeds. The motion of the system
is not calculated by the hydrodynamic system, which is only
concerned with the calculation of the forces.

So to define the hydrodynamic system for the RM3 example, the first
step is to create two hydroBody objects:

.. literalinclude:: /examples/+example_rm3/run_wecsim.m
   :start-at: %% Hydrodynamic body system
   :end-before: % set up transient simulation

The hydroBody objects take as input a file name which is either .h5
file for the whole hydrodynamic system, or a .mat file specific to
that body. The file is always searched for in the hydroData
subdirectory of the case directory, so only the file name is
required, not the full path.

Having created all the required hydroBody objects, we can put them
into a hydroSystem [#f2]_ and tell it to prepare them for a time domain
simulation. You must always call `initialiseHydrobodies` and
`timeDomainSimSetup` before a simulation.

.. literalinclude:: /examples/+example_rm3/run_wecsim.m
   :start-at: % make a hydrosys object for simulation
   :end-before: %% Multibody dynamics system specification

Our hydrodynamic system is now defined, and ready to be simulated,
the next step is then to create the overall multibody dynamics
system which constrains the motion of the bodies in the appropriate
ways. The hydroSystem also takes as input the waveSettings and
simSettings objects, from which it knows the location of the case
directory etc.

.. [#f2] In Matlab we can put the hydroBody objects into the
   hydroSystem simply by placing them in square brackets to make an
   array of hydroBody objects (which is what hydroSystem expects),
   i.e. ``[float_hbody, spar_hbody]``. However, at the time of
   writing, Octave does not yet support this syntax and some small
   modifications are necessary. It simply requires that you instead do::

      obj_array(1) = float_hbody;
      obj_array(2) = spar_hbody;
      hsys = wsim.hydroSystem (waves, simu, obj_array);

   See the file `+example_rm3/run_wecsim_in_octave.m` for an example of
   how the run the same system in Octave. The changes are small, but
   important. This file also demonstrates the use of the .mat file
   hydroData format. This version of the file also runs in Matlab,
   so it is possible to create files which run under both platforms.


Creating the MultiBody System
-----------------------------

Having set up the hydrodynamic system, the next step is to define
the constraints that determine the motion of that system under the
applied forces and moments. These constrained multibody system
dynamics are solved by `MBDyn`_. MBDyn is a separate program, and
forces and motion are communicated between the two programs
throughout the simulation via one of several available
communication methods.

.. _MBDyn: https://www.mbdyn.org/

MBdyn normally takes as input a special data file format which
defines all of the bodies and joints etc. in the system. To make the
creation of systems simpler and allow systems to be completely
defined using Matlab code, an advanced pre-processing tool has been
developed as a part of a suite of tools in the MBDyn Matlab Toolbox,
which is available as a separate standalone tool (and most likely
you will have received it alongside |TNShort|). The MBDyn Matlab
Toolbox also has its own manual with simple examples demonstrating
the toolbox's capabilities. It is advisable to first examine some of
the examples in this manual to get a better understanding of how the
pre-processor works before attempting to develop a system for
|TNShort|, or even understand the example which is about to be
described, as an in-depth discussion of the operation and
organisation of MBDyn and the MBDyn preprocessor will not be
repeated here.

The first step in creating the multibody system is to get the
multibody elements corresponding to the hydrodynamic bodies. This
task can be done automatically by `wsim.hydroSystem` using the
`makeMBDynComponents` method.

.. literalinclude:: /examples/+example_rm3/run_wecsim.m
   :start-at: %% Multibody dynamics system specification (mbdyn)
   :end-at: hsys.makeMBDynComponents ();

This method returns three cell arrays containing the structural
nodes, bodies and other elements which define the parts of the
system corresponding to the hydrodynamic bodies. These components
can then be linked with other components in the system. As the
system requires more than a few lines of code to specify, it has
been specified in a separate (heavily commented) function file in
`+example_rm3/make_multibody_system.m`. A section of the file is
shown below:

.. literalinclude:: /examples/+example_rm3/make_multibody_system.m
   :start-at: % create a node to clamp
   :end-at: % create an orientation with axis 3 pointing along the global axis 2

The full system requires the creation of five joint elements and
another node to act as a reference so the motion can be restricted
to being only in the X-Z plane. This function assembles all of the
multibody dynamics elements into a ``mbdyn.pre.system`` object which
organises and controls them in a similar way to the
``wsim.hydroSystem`` class.

.. literalinclude:: /examples/+example_rm3/make_multibody_system.m
   :start-at: % assemble the system
   :end-before: initptodpos

In the run_wecsim script, it is this system with is returned by the
call to ``example_rm3.make_multibody_system``. The resulting system
can be plotted in a figure:

.. literalinclude:: /examples/+example_rm3/make_multibody_system.m
   :start-at: % draw it in a figure
   :end-before: %% Set up Power Take-Off

The result of the call to ``draw`` above is shown in
:numref:`example_rm3_mbsys_draw_1`

.. _example_rm3_mbsys_draw_1:
.. figure:: /images/example_rm3_mbsys_draw_1.svg
   :align: center

Many options are possible with the draw option, for example, a
prettier output (but less useful for checking the locations of nodes
etc.) can be produced with the 'solid' mode::

   mbsys.draw ( 'Mode', 'solid', ...
                'Light', true, ...
                'AxLims', [-30, 30; -30, 30; -35, 35], ...
                'Joints', false, ...
                'StructuralNodes', false);

The result of this is shown in :numref:`example_rm3_mbsys_draw_2`

.. _example_rm3_mbsys_draw_2:
.. figure:: /images/example_rm3_mbsys_draw_2.svg
   :align: center


Adding a Power Take-Off
-----------------------

The next step is to add a power take-off to the simulation. For this
example we want to add simple damper, i.e. a force which is linearly
related to the relative velocity of the two parts of the RM3 device.
We will calculate this force in Matlab, base on the motion of the
WEC, and send it to MBDyn during the simulation. To implement this,
it is necessary to first get the motion in 3D space of float and
spar and determine their relative velocity in a direction parallel
to the orientation of the spar. Then the damping force must be
calculated, and converted to a 3D force vector to be applied  to the
nodes attached to each body in MBDyn.

|TNShort| makes the process of calculating and applying power
take-off forces easy by providing a set of power take-off classes
(for different types of motion) which perform most of this work,
leaving only the calculation of scalar force value to the user.
These classes are the ``wsim.linearPowerTakeOff`` and
``wsim.rotaryPowerTakeOff`` classes, and their main function is to
automate getting the correct nodal motion from MBDyn during the
simulation, and from this motion determine the relative velocity and
position of two components in a simulation (or the relative angular
velocity and position in the case of ``wsim.rotaryPowerTakeOff``),
and then calculate the force vector on the nodes. Both of these
classes come with extensive help which can be accessed using the
normal Matlab help systems, e.g. run

::

   doc wsim.linearPowerTakeOff

to open the help for this class in the help browser, or

::

   help wsim.linearPowerTakeOff

to view text help in the command line. See
:ref:`required-knowledge-help-system` for more information on the
Matlab and Octave help systems.

For this RM3 example we will use the `wsim.linearPowerTakeOff`
class. This class allows you to apply a force which is calculated by
a function with the syntax:

::

   force = force_function (time, displacement, velocity)

You are free to create and use any Matlab function with this syntax
as the calculation of the power take-off force. For the RM3 example,
the supplied function is a spring-damper function, although in this
case, the spring constant is set to zero. In the example, the
function is defined in the body of the script as an `anonymous
function`_, but it can also be any other function on the Matlab
path. In this case the function accepts the time variable as an
input, but just ignores it.

.. _anonymous function: https://uk.mathworks.com/help/matlab/matlab_prog/anonymous-functions.html

.. literalinclude:: /examples/+example_rm3/run_wecsim.m
   :start-at: %% Set up Power Take-Off
   :end-before: %% Run the simulation

It can be seen that setting up of the PTO is made very simple by
``wsim.linearPowerTakeOff``. It simply requires the two nodes attached
to each part of the PTO, an axis number, and the force function. The
PTO class then does all the work of applying the correct forces to
the nodes. The PTO calculates the forces based on the relative
displacement and velocity of the two nodes in a direction parallel
to the specified axis number in the reference frame of the first node.


Running the Simulation
----------------------

With all of the possible components describing the system now
assembled it would now be possible to start manually stepping
through the system, and indeed this is possible (in case total
control is needed). However, to make running and managing the
simulation easy, the `wsim.wecSim` class has been provided. This
class takes in all the previously define components and runs the
simulation. One of the most useful aspects of `wsim.wecSim` is that
is ensures the correct hydrodynamic and PTO force are applied to the
correct nodes. The `wsim.wecSim` class also does extensive logging
of simulation data for examination afterwards. This logging is also
highly configurable so it can be selectively deactivated to improve
simulation speed. [#fsimlogspeed]_. The logging during a simulation
is performed by another class, the `wsim.logger` class. A logger
object is created inside `wsim.wecSim` and returned at the end of a
simulation. The settings controlling the logging, are themselves
contained in *another* class, `wsim.loggingSettings`. It is the
`loggingSettings` class which is created next in the example script.

.. [#fsimlogspeed] Logging of data usually requires some reallocation
   of memory. ``The wsim.logger`` class used in |TNShort| attempts
   to do this in an efficient way (one can preallocate the size of
   the memory required to log the data, and memory is automatically
   grown in blocks when the preallocated limit is reached). However,
   the time spent doing this is still non-negligible. If efficiency
   is crucial, e.g. within an optimisation procedure, one should
   only log what is really necessary for optimum simulation speed.

.. literalinclude:: /examples/+example_rm3/run_wecsim.m
   :start-at: %% Run the simulation
   :end-before: % create the wecSim object

It can be seen that there is quite fine grained control over what is
logged. It is not necessary to set every setting as shown in the
script, this is just provided to demonstrate what settings are
available. By default, all logging is active (everything will be
logged). One can then selectively deactivate components by setting
them to `false`. Alternatively, the `loggingSettings` class has a
method `allOff`, which deactivates all logging (sets all logging
settings to `false`). One can then selectively reactivate only those
variables which are to be logged during the simulation by setting
them to `true` [#floggingindirect]_.

Having chosen the desired logging settings, the simulation is ready
to be run. The next step is then to create the `wsim.wecSim` object
discussed previously.

.. literalinclude:: /examples/+example_rm3/run_wecsim.m
   :start-at: % create the wecSim object
   :end-before: % initialise the simulation

The `wsim.wecSim` object must take as input the `wsim.hydroSystem`
and `mbdyn.pre.system` objects on creation. The supply of one or
more PTO objects is optional (and they can also be added later using
the `addPTO` method), as is the supply of a `wsim.loggingSettings`
object. If no `loggingSettings` onject is supplied, one is created
internally by `wsim.wecSim` with all default settings, so everything
will be logged. See the help for `wsim.wecSim` for more detail on
how to create the object and what options are available.

Having created the `wecSim` management object we can then prepare
and run the simulation. The `prepare` method must always be called
before a simulation can be run.

.. literalinclude:: /examples/+example_rm3/run_wecsim.m
   :start-at: % initialise the simulation
   :end-before: %% Plot some results

The `run` method then does several things. It creates an input file
for MBDyn and writes it to disk [#foutpuloc]_, it then launches the
MBDyn program and sets up communication with it, and it then steps
though the simulation and returns two outputs. The first output is
the `wsim.logger` object discussed previously which contains the
logged data from the simulation. The second output is an object of
the class `mbdyn.postproc`. Both of these outputs will be discussed
in more detail in the next section.

There are many options available for the `run` method which are
documented in it's help. In this case we used only one, the
``'TimeExecution'`` option, which prints the wall-clock time taken
to perform the simulation.

.. [#floggingindirect] It should also be noted that setting some
   logging settings to `true` indirectly causes some other variables
   to also be logged. For example, setting `forceAddedMass` to
   `true` will also cause `forceAddedMassUncorrected` to become
   `true` as this is required to calculate the final added mass force.

.. [#foutpuloc] By default these files are written to a subdirectory
   of the case directory named `output_<date_and_time>` where
   `<date_and_time>` is replaced by the date and time the run method
   was started. The output from MBDyn is also placed in this directory.


Examining The Results
---------------------

Having run the simulation, it is now possible to use the returned
`wsim.logger` object and `mbdyn.postproc` object to examine what
occurred during the simulation. The `wsim.logger` object generally
contains data which was created on the Matlab side of the
simulation, i.e. the hydrodynamic and PTO forces, but also data
which was sent by MBDyn during the simulation such as node positions
etc. The `mbdyn.postproc` method in contrast can be used to access
other more detailed data calculated by MBDyn during the simulation
and written to disk once the simulation was complete. This includes
data on all the internally calculated forces etc.

The data in `wsim.logger` object is stored in a public property,
`data`, which is a Matlab structure with field names corresponding
to each logged data item. For example, the contents of `data` in the
`datalog` output of the RM3 example is the following:

.. highlight:: none

::

   >> datalog.data

   ans =

     struct with fields:

                             Time: [4001×1 double]
                        Positions: [3×2×4001 double]
                 AngularPositions: [3×2×4001 double]
                       Velocities: [3×2×4001 double]
                AngularVelocities: [3×2×4001 double]
                    Accelerations: [3×2×4001 double]
             AngularAccelerations: [3×2×4001 double]
                       NodeForces: [3×2×4001 double]
            NodeForcesUncorrected: [3×2×4001 double]
                       ForceHydro: [3×2×4001 double]
                  ForceExcitation: [3×2×4001 double]
              ForceExcitationRamp: [3×2×4001 double]
               ForceExcitationLin: [3×2×4001 double]
            ForceExcitationNonLin: [3×2×4001 double]
            ForceRadiationDamping: [3×2×4001 double]
                   ForceRestoring: [3×2×4001 double]
                    ForceMorrison: [3×2×4001 double]
              ForceViscousDamping: [3×2×4001 double]
                   ForceAddedMass: [3×2×4001 double]
        ForceAddedMassUncorrected: [3×2×4001 double]
                      NodeMoments: [3×2×4001 double]
           NodeMomentsUncorrected: [3×2×4001 double]
                      MomentHydro: [3×2×4001 double]
                 MomentExcitation: [3×2×4001 double]
             MomentExcitationRamp: [3×2×4001 double]
              MomentExcitationLin: [3×2×4001 double]
           MomentExcitationNonLin: [3×2×4001 double]
           MomentRadiationDamping: [3×2×4001 double]
                  MomentRestoring: [3×2×4001 double]
                   MomentMorrison: [3×2×4001 double]
             MomentViscousDamping: [3×2×4001 double]
                  MomentAddedMass: [3×2×4001 double]
       MomentAddedMassUncorrected: [3×2×4001 double]
              PTO_1_InternalForce: [4001×1 double]
       PTO_1_RelativeDisplacement: [4001×1 double]
           PTO_1_RelativeVelocity: [4001×1 double]

.. highlight:: default

This can be accessed an processed or plotted like any normal Matlab
variable. However, for convenience, the `wsim.logger` class provides
methods to plot the outputs directly. The main method for this is
`plotVar`. The `plotVar` method is called with the name of the
variable to be plotted and plots it against it's independent
variable, which in most cases is the *Time* variable. This is shown
for several variables in the example:

.. literalinclude:: /examples/+example_rm3/run_wecsim.m
   :start-at: %% Plot some results
   :end-before: % we can also plot the motion

The plots produced by these calls to `plotVar` are shown in
:numref:`rm3_positions`, :numref:`rm3_velocities` and
:numref:`rm3_pto_force`.

.. _rm3_positions:
.. figure:: /images/rm3_positions_output.png
   :align: center

.. _rm3_velocities:
.. figure:: /images/rm3_velocities_output.png
   :align: center

.. _rm3_pto_force:
.. figure:: /images/rm3_pto_force_output.png
   :align: center

The log returned by the `wecSim` object only contains motion data
for the nodes which are accessed through the external structural
force, there may be other nodes in the simulation, and data on the
other nodes in the simulation must be obtained from the output of
MBDyn. As mentioned previously, the second output of the `run`
method, is the `mbdyn.postproc` object. This object loads data from
a netcdf format data file produced by MBDyn at the end of the
simulation. A list of all the variables which can then be examined
using this object is produced using the `displayNetCDFVarNames`
method, e.g. for the RM3 example, this produces the following::

   >> mbdyn_pproc.displayNetCDFVarNames ()
   run.step : time step index
   time : simulation time
   run.timestep : integration time step
   node.struct : Structural nodes labels
   node.struct.1 : no description
   node.struct.1.X : global position vector (X, Y, Z)
   node.struct.1.R : global orientation matrix (R11, R21, R31, R12, R22, R32, R13, R23, R33)
   node.struct.1.XP : global velocity vector (v_X, v_Y, v_Z)
   node.struct.1.Omega : global angular velocity vector (omega_X, omega_Y, omega_Z)
   node.struct.2 : no description
   node.struct.2.X : global position vector (X, Y, Z)
   node.struct.2.R : global orientation matrix (R11, R21, R31, R12, R22, R32, R13, R23, R33)
   node.struct.2.XP : global velocity vector (v_X, v_Y, v_Z)
   node.struct.2.Omega : global angular velocity vector (omega_X, omega_Y, omega_Z)
   node.struct.3 : no description
   node.struct.3.X : global position vector (X, Y, Z)
   node.struct.3.R : global orientation matrix (R11, R21, R31, R12, R22, R32, R13, R23, R33)
   node.struct.3.XP : global velocity vector (v_X, v_Y, v_Z)
   node.struct.3.Omega : global angular velocity vector (omega_X, omega_Y, omega_Z)
   elem.autostruct : AutomaticStructural elements labels
   elem.joint : Joint elements labels
   elem.force : Force elements labels
   node.struct.1.B : momentum (X, Y, Z)
   node.struct.1.G : momenta moment (X, Y, Z)
   node.struct.1.BP : momentum derivative (X, Y, Z)
   node.struct.1.GP : momenta moment derivative (X, Y, Z)
   node.struct.2.B : momentum (X, Y, Z)
   node.struct.2.G : momenta moment (X, Y, Z)
   node.struct.2.BP : momentum derivative (X, Y, Z)
   node.struct.2.GP : momenta moment derivative (X, Y, Z)
   elem.joint.7 : no description
   elem.joint.7.f : local reaction force (Fx, Fy, Fz)
   elem.joint.7.m : local reaction moment (Mx, My, Mz)
   elem.joint.7.F : global reaction force (FX, FY, FZ)
   elem.joint.7.M : global reaction moment (MX, MY, MZ)
   elem.joint.8 : no description
   elem.joint.8.f : local reaction force (Fx, Fy, Fz)
   elem.joint.8.m : local reaction moment (Mx, My, Mz)
   elem.joint.8.F : global reaction force (FX, FY, FZ)
   elem.joint.8.M : global reaction moment (MX, MY, MZ)
   elem.joint.10 : no description
   elem.joint.10.f : local reaction force (Fx, Fy, Fz)
   elem.joint.10.m : local reaction moment (Mx, My, Mz)
   elem.joint.10.F : global reaction force (FX, FY, FZ)
   elem.joint.10.M : global reaction moment (MX, MY, MZ)
   elem.joint.10.R : global orientation matrix (R11, R21, R31, R12, R22, R32, R13, R23, R33)
   elem.joint.10.Omega : local relative angular velocity (x, y, z)


This shows the list of available variables and a short description
of each. The contents of these can then be accessed using the
`getNetCDFVariable` method. When it is created, the `mbdyn.postproc`
class also immediately loads the motion data from the MBDyn output
file for all the nodes in the simulation. Some methods for plotting
this motion data are then provided, such as the
`plotNodeTrajectories` method demonstrated in the script:

.. literalinclude:: /examples/+example_rm3/run_wecsim.m
   :start-at: %% Plot some results
   :end-before: % we can also plot the motion

The output of this command is shown in :numref:`rm3_postproc_traj`.
As usual you can learn about the full range of post-processing
methods available from the `mbdyn.postproc` object by examining the
help for the class, but also through the documentation for the
Matlab NBDyn Toolbox.

.. _rm3_postproc_traj:
.. figure:: /images/postproc_node_trajectories.png
   :align: center

Finally, it can also be extremely helpful to visualise the motion of
the system during the simulation. This can be done using the
`animate` method of the `wsim.wecSim` class. An example of this is
shown at the end of the script.

.. literalinclude:: /examples/+example_rm3/run_wecsim.m
   :start-at: %% Animate the system

A still image from the resulting animation is shown in
:numref:`rm3_animate`. It can be seen that the wave surface is
plotted in addition to the device. This animation by defuat is just
played in a Matlab figure window, however, if desired, it can also
be written to disk as an avi file. See the help for the
`wsim.wecSim.animate` for more information.

.. _rm3_animate:
.. figure:: /images/rm3_animation_still.png
   :align: center

Conclusions
===========

This example gives a good introduction to the capabilities of
|TNShort|, but is not exhaustive. Topics which of interest which
have not been explored include multi-rate simulation techniques,
creating your own power take-off classes derived from the built-in
power take-off classes, and advanced features of MBDyn that may be
useful for system simulation.

Note also that further examples are provided in the same location as
the example described in this document.
