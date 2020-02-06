.. Renewnet Foundry documentation master file, created by
   sphinx-quickstart on Fri May  5 14:45:20 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Renewnet Foundry
****************

The Renewnet Foundry is a suite of tools originally created for the
simulation of renewable energy devices and systems. It is developed
in Matlab, but most componenets can also work in Octave. However it
is also useful for simulating many other kinds of systems and
subsystems, as it is a modular suite of tools with a range of
capabilities. The most well developed parts of the foundry are the
following tools:

* MBDyn Matlab Toolbox
* Edinburgh Wave Systems Toolbox
* Edinburgh Electrical Machines Toolbox

In addition to these main tools, there are also Matlab functions and
classes, including plotting tools, utilites for checking inputs,
exporting data and much more.

INSTALLATION
============

First download the release .zip file and extract it's contents on
your computer.

.. note:: Windows users, make sure you follow this step, sometimes it is easy to forget you are just browsing the still-compressed zip folder in the file browser.

To get started with the matlab code you should run the function::

  rnfoundry_setup

which you can find in the top level directory. It is also advised that you
READ THE HELP FOR THIS FUNCTION before running it to get an idea
of what it does and what the system requirements are to get
optimum performance.

DOCUMENTATION
=============

Documentation is provided in two forms, inside the matlab functions
and classes, in the standard Matlab help format, and for some tools,
as standalone documentation in html (web page) format. The html
documentation can be found in the *documentation* subdirectory in the
top level directory. The latest version of the documentation for some
of the tools is also hosted online, e.g:

* |ewst-link-pre|\ |version|\ |ewst-link-mid|\ |version|\ |ewst-link-post|
* |mbdyn-link-pre|\ |version|\ |mbdyn-link-mid|\ |version|\ |mbdyn-link-post|

.. |ewst-link-pre| raw:: html

    <a href="documentation/EWST-

.. |ewst-link-mid| raw:: html

    -html-docs/index.html">Edinburgh Wave Systems Toolbox Documentation v

.. |ewst-link-post| raw:: html

    </a>

.. |mbdyn-link-pre| raw:: html

    <a href="mbdyn-matlab-toolbox-

.. |mbdyn-link-mid| raw:: html

    -html-docs/index.html">MBDyn Matlab Toolbox Documentation v

.. |mbdyn-link-post| raw:: html

    </a>

USAGE
=====

You can find some example scripts in the following directories:

Edinburgh Wave Systems Toolbox
------------------------------

rnfoundry/wave/doc/sphinx/examples

This contains examples of using the Edinburgh Wave Systems Toolbox.
Further (less refined and well commented) examples of using this may
be found in

rnfoundry/wave/matlab-octave/wec-sim/test


MBdyn Multibody Dynamics Toolbox
--------------------------------

rnfoundry/common/multibody/MBDyn/doc/sphinx/examples

Further (less refined and commented) examples can be found in the
testing code in

rnfoundry/common/multibody/MBDyn/test


Permanent Magnet Machines Toolbox
---------------------------------

rnfoundry/common/electrical/matlab-octave/permanent_magnet_machines_tools/examples_and_tutorials

Contains examples of using the permanent magnet machine simulation
tools


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   self



Indices and tables
******************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
