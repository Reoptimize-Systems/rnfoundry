Installation
************

All |TNShort| functionality is based on Matlab code, as such installation 
requires only that the code is present on the `Matlab path`_. 
However, it should be noted that |TNShort| is dependant on the MBDyn 
Matlab preprocessor, which has its own installation instructions and 
requirements. It is most likely you will have received both EWST and 
the required MBDyn tools as a single package. These two packages, 
while standalone, are closely related and developed in tandem.

.. _Matlab path: https://uk.mathworks.com/help/matlab/matlab_env/what-is-the-matlab-search-path.html


RenewNet Foundry
================

You may also have received |TNShort| as part of the RenewNet Foundry 
package, in this case you should first refer to the README file in 
the top level of the RenewNet Foundry package and follow the 
instructions there. The RenewNet Foundry installation is automated 
through the `rnfoundry_setup.m` function, and |TNShort| is installed 
automatically by this.

.. _install-nemoh:

Nemoh
=====

The toolbox provides an interface to the Nemoh hydrodynamics BEM 
solver used for calculating wave loads on offshore structures. A 
copy of the Nemoh solver is not provided with the toolbox. You must 
obtain this software separately. Further information on where Nemoh 
should be installed is available in the section 
_nemoh-interface-exec-loc. For convenience, a copy of Nemoh 2.03 
(the latest stable release) may be downloaded from the following 
location: `RenewNet Foundry Nemoh Downloads`_. 

.. _RenewNet Foundry Nemoh Downloads : https://sourceforge.net/projects/rnfoundry/files/Nemoh/

Instruction and information about the Nemoh download are provided on 
this page.
