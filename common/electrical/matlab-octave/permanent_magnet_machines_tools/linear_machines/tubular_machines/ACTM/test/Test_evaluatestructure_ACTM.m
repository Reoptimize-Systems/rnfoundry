% Test_evaluateStructure_ACTM

Rm = 0.05;
g = 6/1000;
Ri = Rm + g;
WmVWp = 0.5;
WpVRm = 0.5;
RiVRm = Ri / Rm;
RoVRm = 1.1;
RaVRo = 1.1;
RsoVRm = 0.4;
cwVWp = 1/3;
kfill = 0.65;
dc = 1e-4;        % 1 mm diameter wire 
RsiVRm =  0; 
N = 2867;
supportLengths = [2, 2; 2, 2]; 
totalLength = [19.125, 19.125];
sections = 10;
vCoil = 0.3;
E = [207000000000.000,10000000,207000000000.000;];
Sz = 12;
Sr = 8; 
steelDensity = 7800;
magDensity = 7500;
shaftDensity = 7500;
copperDensity = 8960;
% avAxBVals = 


[PorF, z, actForce, maxDef, avAxBVals] = evaluatestructure_ACTM(WmVWp, WpVRm, RoVRm, RaVRo, RsiVRm, RsoVRm,...
                                                                cwVWp, Rm, kfill, g, N, dc, supportLengths,...
                                                                totalLength, sections, vCoil, E, Sz, Sr, steelDensity,...
                                                                magDensity, shaftDensity, copperDensity);