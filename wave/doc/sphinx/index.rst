.. EWST documentation master file, created by
   sphinx-quickstart on Fri May  5 14:45:20 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

|TNShort| (|TNFull|)
****************************************

The |TNShort| (|TNFull|) is a suite of tools for the
simulation of wave energy devices. |TNShort| has the ability to model devices that
are comprised of rigid bodies, power-take-off systems, and mooring systems.
Simulations are performed in the time-domain by solving the governing WEC
equations of motion in 6 degrees-of-freedom.

|TNShort| is derived from `WEC-Sim`_, an open-source wave energy converter simulation tool developed in Matlab/Simulink
using the multi-body dynamics solver Simscape Multibody. The |TNShort| is also
developed in Matlab, but drops the requirement for Simulink or Simscape
Multibody. It also aims to be compatible with `Octave`_,
an alternative system able to process much of the standard Matlab code base.
The hydrodynamic data format for both is identical, so hydrodynamic data can
be easily ported between systems.

|TNShort| replaces the mutlibody modelling parts of the code with 
`MBDyn`_, an advanced multiboy dynamics simulator. An advanced 
Matlab code based MBDyn preprocessor is available to allow the 
creation of MBDyn model input files directly from Matlab. Detailed 
documentation of the preprocesing tool may be found in it's own 
dedicated manual.

.. _WEC-Sim: http://wec-sim.github.io/WEC-Sim/index.html
.. _Octave: https://www.gnu.org/software/octave/
.. _MBDyn: https://www.mbdyn.org/

Purpose
=======

With the existance of Wec-Sim the need for |TNShort| may nt be obvious. |TNShort| addresses
several main perceived issues with Wec-Sim.

1. *Maintainability:* Wec-Sim being heavily based in Simulink does not lend itself
to maintenance using standard software version control systems (e.g. `Git`_,
`Mercurial`_). Simulink files cannot be easily compared
for changes using these systems.

2. *Debugability:* Debugging a purely code based system is easier than one based
in Simulink, as the ability to step through the code and navigate the different levels
of the system is more advanced.
purely code-based

3. *Interface design:* The developers of |TNShort| have a different interface design
philosophy which is more conventional than the WEC-Sim method. The |TNShort| interface
is more oriented toward automation and batch processing than the original WEC-Sim
interface. This is mainly to facilitate randomised simulation and optimisation
algorithms.

4. *Modifiability:* Being purely code based, and with all code in one location
(rather than spread throughout Simulink models)

5. *Cost:* WEC-Sim requires many commercial platforms to be operated, |TNShort| can be run
entirely on free software (although the performance on Matlab will be
superior).

The |TNShort| developers recognise that not everyone will agree with the points above
or that they justify the creation of a separate system, but it was these needs which
drove it's creation.

An aim of the project will continue to be to maintain compatibility as much as
possible between the two systems

.. _Git: https://git-scm.com/
.. _Mercurial: https://www.mercurial-scm.org/

Developers and Contributers
===========================

The |TNFull| has been created by The Institute for Energy Systems at The University
of Edinburgh as part of the EPSRC funded project "All Electrical Drive Train for
Marine Energy Converters (EDRIVE-MEC)", grant No. EP/N021452/1. The main contributor
is Dr. Richard Crozier.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   self
   installation
   getting_started
   required_knowledge
   nemoh
   api_reference



Indices and tables
******************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
