% Test_AirGapClosure_ACPMSM: a script to test the function
% airgapclosure_ACPMSM
%
%
% gapclosingforce_ACPMSM(1.93, 8.708, 1.3, 4.4, 1.5, 0.09)
E = 200e9;


dgVlm = 2.0;
bpVlm = 8.7;
taupVbp = 1.3;
lsVbp = 4.4;
dbiVlm = 1.5;
lm = 0.02;
sections = 10;

%           IVars(1,1): b, width of the I-beam flanges
%           IVars(1,2): t, the thickness of the flanges
%           IVars(1,3): tw, the thickness of the I-Beam vertical part (Web)
%           IVars(1,4): d, height of the I-Beam vertical part
bp = bpVlm * lm;

IVars = [(taupVbp*bp) (30/1000) (8.4/1000) (426/1000)];

dgThresh = 3/1000;

dg = airgapclosure_ACPMSM(E, IVars, dgVlm, bpVlm, taupVbp, lsVbp, dbiVlm, lm, sections, dgThresh)