Getting Started
***************

The MBDyn Toolbox can be used to generate input files for MBDdyn, 
run the simulation, and load the results into Matlab. In addition, 
co-simulation, where forces are sent to MBDyn from Matlab and motion
data sent to Matlab from MBDyn on each time step is also possible.

The MBDyn solver takes as input a text-based description of the 
multibody dynamics problem. The format of the MBDyn input files can

The MBDyn pre-processing tools which generate the files are based on
an object oriented appoach, elements of the problem are represented as 
classes and a special system class acts as a container for the entire 
collection of input data, conceptually this is shown in 
:numref:`mbdyn_basic_object_diagram`.

.. _mbdyn_basic_object_diagram:
.. figure:: /images/mbdyn_basic_object_diagram.svg



