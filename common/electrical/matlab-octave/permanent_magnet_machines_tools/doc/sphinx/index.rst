.. EMST documentation master file, created by
   sphinx-quickstart on Fri May  5 14:45:20 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. highlight:: matlab

EMST (Electrical Machine Simulation Toolchain)
**********************************************

EMST (Electrical Machine Simulation Toolchain) is a suite of tools for the
simulation of electrical machines, particularly permanent magnet machines.

EMST supports the simulation and analysis of a wide variety of machine types,
including linear and rotary machine types, and for each of these, air-cored,
slotless and slotted variations. For rotary machines both radial flux and axial
flux machine types are available, and for linear machines, both conventional
flat profiled machines and tubular configurations are possible.

In addition to this, the machine models are constructed in a modular way
with code that can be reused to easily add new machine types from more basic
components, e.g. using an existing armature drawing code with a new type of
permanent magnet rotor not currently available, and also reusing the many
utility functions such as winding layout calculations, iron loss
calculation codes, dynamic simulation tools etc.

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
