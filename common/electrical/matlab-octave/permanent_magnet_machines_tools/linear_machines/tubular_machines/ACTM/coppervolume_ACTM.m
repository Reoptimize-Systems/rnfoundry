function CopperVol = coppervolume_ACTM(RoVRm, Ntot, Rm, g, dc, armPoles)
% coppervolume_ACTM: a function for calculating the total volume of copper
% in the ACTM armature. Round wire is asumed
%
% Input: 
%
%   RmVRo - Translator radius to coil outer radius ratio
%
%   Ntot - number of turns in a single coil
%
%   Rm - translator radius in m
%
%   g - air-gap in m
%
%   dc - wire diameter in m
%
%   armPoles - number of Poles in the armature
%
% Output:
%
%   CopperVol - total volume of copper used in the machine in m^3

    CopperVol = 3 * armPoles * wirelength_TM(RoVRm, Ntot, Rm, g) * pi * (dc/2)^2;

end