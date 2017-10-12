.. MBDyn interface documentation master file

MBDyn Matlab Toolbox
********************

The MBDyn Matalb Toolbox is a suite of tools for interacting with 
MBDyn from matlab.

`MBDyn`_ is a free Multibody Dynamics analysis software, released under 
GNU's GPL 2.1 developed at the Dipartimento di Scienze e Tecnologie 
Aerospaziali (formerly Dipartimento di Ingegneria Aerospaziale) of the 
University "Politecnico di Milano", Italy.

.. _mbdyn: https://www.mbdyn.org/

MBDyn features the integrated multidisciplinary simulation of 
multibody, multiphysics systems, including nonlinear mechanics of 
rigid and flexible bodies (geometrically exact & composite-ready 
beam and shell finite elements, component mode synthesis elements, 
lumped elements) subjected to kinematic constraints, along with 
smart materials, electric networks, active control, hydraulic 
networks, and essential fixed-wing and rotorcraft aerodynamics.

MBDYn does not feature a GUI and requires that problems are 
described using a text based input files. The MBDyn Matlab Toolbox 
includes a preprocessor for generating these files using Matlab code, 
with visualisation tools to assist development, an interface for 
communication of forces and motions between Matlab and MBDyn during 
a simulation and a post-processor with visualisation and animation 
capabilities.

Octave
======

The MBDyn Matlab Toolbox aims to also target `Octave`_, and most parts
of the toolbox will also be supported on this platform. However, 
there are some limitations. The main area where Octave is not fully 
compatible is the visualisation tools. Some of these tools currently use
functions and classes not available in Octave at the time of writing. 
The Octave developers aim for compatibility with Matlab, so it is expected  
that the functions will become available in due course. Otherwise,  
throughout this documentation you may consider "Matlab" to be 
interchangeable with "Octave".

.. _Octave: https://www.gnu.org/software/octave/


Developers and Contributers
===========================

The MBDyn Matlab Toolbox has been created by The Institute for 
Energy Systems (IES) at The University of Edinburgh as part of the 
EPSRC funded project "All Electrical Drive Train for Marine Energy 
Converters (EDRIVE-MEC)", grant No. EP/N021452/1. The main contributor
is Dr. Richard Crozier.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   self
   installation
   getting_started


Indices and tables
******************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
