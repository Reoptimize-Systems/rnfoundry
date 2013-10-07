%% Test_beamdef_PMSM
%
% A script for testing the beamdef_PMSMS function
%
Taup = 0.12;
ls = 1; 
Phases = 3;
Taus=Taup/Phases;
bs=Taus/2; 
bt=Taus/2;
hs=5*bt;

%   IVars - (2 x 4) matrix of values for calculating the second moment of
%           area of the stator and translator sections, the second row of
%           values is used for both stator sections if a double sided
%           machine is being investigated.
%
%           IVars(1,1): b, width of the I-beam flanges
%           IVars(1,2): t, the thickness of the flanges
%           IVars(1,3): tw, the thickness of the I-Beam vertical part
%           IVars(1,4): d, height of the I-Beam vertical part
%
%           IVars(2,1): bc, thickness of central section
%           IVars(2,2): Taup, the pole width
%           IVars(2,3): dt, the tooth height
%           IVars(2,4): bt, the tooth width
IVars = [Taup (Taup/10) (Taup/15) (Taup/2); 0.02 Taup hs bt];

%% Linearly distributed forces
FVars = [-1000 -1000; 1000 -1000; 1000 1000];

totalLength = ls;
aL = ls;
x = 0:totalLength/10:totalLength;
E = [207e9 151e9];

beamdef_PMSM(IVars, totalLength, aL, x, FVars, E);

%% non-linearly Distributed forces
FVars = [-100 -200 -560 -2000 -2000 -560 -200 -100;...
          0    0    0    0     0     0    0    0;...
          100  200  560  2000  2000  560  200  100];
      
%aL = [0  0.125 0.25 0.375  0.5  0.625 0.75 0.875 1];
aL = 0:totalLength/size(FVars,2):totalLength;

totalLength = ls;

x = 0:totalLength/10:totalLength;
E = [207e9 151e9];

Def = beamdef_PMSM(IVars, totalLength, aL, x, FVars(1:2,:), E);
