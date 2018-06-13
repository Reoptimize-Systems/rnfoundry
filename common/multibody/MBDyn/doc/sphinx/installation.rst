Installation
************

The MBDyn Toolbox consists of normal Matlab code files, and one C++ mex 
interface file which must be compiled before it can be used in Matlab. 
Releases will be distributed with the compiled mex file for several 
platforms, so if yusing an official release, you will not need to do 
any compilation and all that is required is that the code be present on 
the `Matlab path`_,  

.. _Matlab path: https://uk.mathworks.com/help/matlab/matlab_env/what-is-the-matlab-search-path.html

You can add the code to the Matlab path by using the "Set Path" button 
which opens the Set Path dialog. From here, choose "Add With Subfolders"
and select the MBDyn code directory. 

Alternatively you can add the paths using the Matlab ``addpath`` and 
``genpath`` commands, e.g.::

    addpath ( genpath ('/my/code/directory/MBDyn') );

See the help for ``addpath`` and ``genpath`` for more information.

RenewNet Foundry
================

You may also have received the MBDyn Toolbox as part of the RenewNet 
Foundry package, in this case you should first refer to the README 
file in the top level of the RenewNet Foundry package and follow the 
instructions there. The RenewNet Foundry installation is automated 
through the `rnfoundry_setup.m` function, and the MBDyn Toolbox is 
installed automatically by this.

MBDyn
=====

The MBDyn toolbox naturally requires a copy of MBDyn. Releases of 
the toolbox include a copy of MBDyn, but the official version may 
also be used. However, building MBDyn for Windows is not trivial, 
and official Windows releases are not available, so the bundled 
version is recommended for Windows. The MBDyn toolbox will locate 
the bundled version without help. On linux or unix systems you 
should follow the normal installation instructions for MBDyn. You 
can obtain the latest version of MBDyn from `the MBDyn website`_.

.. _the MBDyn website: https://www.mbdyn.org/?Software_Download

Compiling the mex interface file
================================

If you do not have the mex file for your computer platform, or you want 
to use the development code rather than a release, you will need to 
build it from source. However, the mex interface is stable and unlikely 
to change, so the release version of the compiled mex file is likely to 
work with newer Matlab sources. 

To build the mex file you will also need a `Matlab compatible C++ compiler`_ 
and `set it up within Matlab for use as a mex compiler`_. Development
generally occurs on Linux, and only the gcc compiler on Linux, and 
Mingw compiler on Windows (which is a Windows port of gcc) have been 
tested. The use of Mingw if using Windows is therefore recommended. 

.. _Matlab compatible C++ compiler: https://www.mathworks.com/support/compilers.html
.. _set it up within Matlab for use as a mex compiler: https://uk.mathworks.com/help/matlab/matlab_external/choose-c-or-c-compilers.html

Once a compiler is set up, the mex file can be created it can be 
created by using the function ``mexmbdyn_setup`` found in the top level 
MBDyn directory. See the help for ``mexmbdyn_setup`` for details on the use 
of this function.





