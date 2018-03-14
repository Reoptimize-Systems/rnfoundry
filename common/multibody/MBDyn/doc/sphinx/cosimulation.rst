************
Cosimulation
************

The previous examle simply set up a simulation and ran it, with not 
input from Matlab during the solving process. MBDyn is also capable 
of cosimulation with external software. This allows forces to be 
sent to MBDyn based on kinematics sent by MBDyn to the external 
software. The MBDyn Matlab toolbox includes tools to assist with 
connecting to an MBDyn simulation and sending forces as the 
simulation progresses. Here we will demonstrate the use of these 
tools with a simple example.

The cosimulation process is primarily handled by the 
``mbdyn.mint.MBCNodal`` class which manages the communication with 
MBDyn.
