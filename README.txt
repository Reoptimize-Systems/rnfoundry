RENEWNET FOUNRDY
================

The renewnet foundry is a set of matlab codes and other software
aimed at the simulation of renewable energy systems, particularly
from the perspective of the development of the power take-off
components. The most well developed parts of the foundry are the
permanent magnet machine design and simulation tools.

Most code in the foundry is also compatible with Octave, the free
alternative to Matlab.

INSTALLATION
------------

To get started with the matlab code you should run the function

rnfoundry_setup.m

found in the top level directory. It is also advised that you
READ THE HELP FOR THIS FUNCTION before running it to get an idea
of what it does and what the system requirements are to get
optimum performance. Ideally you will have a C++ compiler set up
with your matlab installation, see the rnfoundry_setup help for
more information. Other than this the setup is fully automated.

USAGE
-----

You can find some example scripts in the following directories:

rnfoundry-hg/wave/matlab-octave/point_absorber/examples

(contains examples of using the heavng buoy simulation)

rnfoundry-hg/common/electrical/matlab-octave/permanent_magnet_machines_tools/examples_and_tutorials

(contains examples of using the permanent magnet machine
simulation tools)
