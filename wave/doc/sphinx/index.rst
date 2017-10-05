.. EMST documentation master file, created by
   sphinx-quickstart on Fri May  5 14:45:20 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. highlight:: matlab

EWST (Edinburgh Wave Simulation Toolbox)
****************************************

The EWST (Edinburgh Wave Simulation Toolbox) is a suite of tools for the
simulation of wave energy devices.

EWST is derived from `WEC-Sim <http://wec-sim.github.io/WEC-Sim/index.html>`, 
an open-source wave energy converter simulation tool developed in Matlab/Simulink
using the multi-body dynamics solver Simscape Multibody. The EWST is also 
developed in Matlab, but drops the requirement for Simulink or Simscape 
Multibody. It also aims to be compatible with `Octave <https://www.gnu.org/software/octave/>`, 
an alternative system able to process much of the standard Matlab code base.
The hydrodynamic data format for both is identical, so hydrodynamic data can
be easily ported between systems. 

EWST replaces the mutliboy modelling parts of the code with the `MBDyn <https://www.mbdyn.org/>`


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   self
   installation
   getting_started
   machine_types
   lower_level_functions




Indices and tables
******************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
