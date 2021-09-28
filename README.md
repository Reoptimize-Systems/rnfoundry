# RENEWNET FOUNRDY

The renewnet foundry is a set of matlab codes and other software
aimed at the simulation of renewable energy systems, including
development of the power take-off components. The most well
developed parts of the foundry are the multibody dynamics tools,
wave systems tools and the permanent magnet machine design and
simulation tools.

Most code in the foundry is also compatible with Octave, the free
alternative to Matlab.

## INSTALLATION

To get started with the matlab code you should run the function::

  rnfoundry_setup

found in the top level directory. It is also advised that you
READ THE HELP FOR THIS FUNCTION before running it to get an idea
of what it does and what the system requirements are to get
optimum performance. Ideally you will have a C++ compiler set up
with your matlab installation, see the rnfoundry_setup help for
more information. Other than this the setup is fully automated.

## USAGE

You can find some example scripts in the following directories:

### Edinburgh Wave Systems Toolbox

rnfoundry/wave/doc/sphinx/examples

This contains examples of using the Edinburgh Wave Systems Toolbox.
Further (less refined and well commented) examples of using this may
be found in

rnfoundry-hg/wave/matlab-octave/wec-sim/test


### MBdyn Multibody Dynamics Toolbox

rnfoundry/common/multibody/MBDyn/doc/sphinx/examples

Further (less refined and commented) examples can be found in the
testing code in

rnfoundry/common/multibody/MBDyn/test


### Permanent Magnet Machines Toolbox

rnfoundry/common/electrical/matlab-octave/permanent_magnet_machines_tools/examples_and_tutorials

Contains examples of using the permanent magnet machine simulation
tools
