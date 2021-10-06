%% Test_BeamDef_ACTIAM
%
% A script for testing the beamdef_ACTIAM function
%

Ra = 0.16;
Ro = 0.15;
Ri = 0.103;
Rm = 0.1;
Rb = 0.0025;
Rso = Rm;
Rsi = Ro - 0.01;

%   IVars - (3 x 2) matrix of values for calculating the second moment of
%           area of the stator and translator sections, the second row of
%           values is used for both stator sections if a double sided
%           machine is being investigated.
%
%           IVars(1,1): Rso, outer radius of shaft support
%           IVars(1,2): Rsi, inner radius of shaft support
% 
%           IVars(2,1): Rm, outer radius of field
%           IVars(2,2): Rb, radius of bore hole in field for tension cable
%
%           IVars(3,1): Ro, outer coil radius
%           IVars(3,2): Ri, inner coil radius
%
%           IVars(4,1): Ra, outer sheath radius
%           IVars(4,2): Ro, inner sheath radius
IVars = [Rso Rsi; Rm Rb; Ro Ri; Ra Ro];

supportLengths = [1 2];
Lf = 3;
La = sum(supportLengths) + Lf;
E = [207e9 151e9 120e9 207e9];
x = 0:La/200:La;

%% Uniform Force

fFVars = 1000;
aFVars = 1000;
Def = beamdef_ACTIAM(IVars, supportLengths, Lf,  La, fFVars, aFVars, E, x);

%% Linearly distributed forces
fFVars = [-1000 1000];
aFVars = [1000 -1000];

Def = beamdef_ACTIAM(IVars, supportLengths, Lf,  La, fFVars, aFVars, E, x);

%% non-linearly Distributed forces
FVars = [-100 -200 -560 -2000 -2000 -560 -200 -100;...
          0    0    0    0     0     0    0    0;...
          100  200  560  2000  2000  560  200  100];
      
%aL = [0  0.125 0.25 0.375  0.5  0.625 0.75 0.875 1];
aL = 0:totalLength/size(FVars,2):totalLength;

totalLength = ls;

x = 0:totalLength/10:totalLength;
E = [207e9 151e9];

Def = BeamDef_PMSM(IVars, totalLength, aL, x, FVars(1:2,:), E);
