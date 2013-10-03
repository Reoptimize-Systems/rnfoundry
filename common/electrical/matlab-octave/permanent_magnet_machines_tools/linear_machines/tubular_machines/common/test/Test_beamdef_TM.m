% Test_BeamDef_TM

%           IVars(1,1): Rso, outer radius of shaft support
%           IVars(1,2): Rsi, inner radius of shaft support
%
%           IVars(2,1): Ro, outer coil radius
%           IVars(2,2): Ri, inner coil radius
%
%           IVars(3,1): Ra, outer sheath radius
%           IVars(3,2): Ro, inner sheath radius
Rso = 0.02;
Rsi = 0.01;
Ro = 0.12;
Ri = 0.1003;
Ra = 0.13;
     
IVars = [Rso Rsi;...
         Ro  Ri;...
         Ra  Ro];
     

supportLengths = [1.0 1.0; 0.5 0.5];

fLength = 3;
aLength = 2;
totalLength = [fLength aLength]; 

aL = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];

P = 10000;

aFVars = [repmat([100],1,size(aL,2)-1)];

fFVars = -aFVars;

E = [200e9 100 200e9];

x = 0:0.05:3;

Def = beamdef_TM(IVars, supportLengths, totalLength, aL, P, fFVars, aFVars, x, E)