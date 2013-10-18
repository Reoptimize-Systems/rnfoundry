% Test_Lorentz_Force_Calc_ACTM

% First set the machine physical variables
design.Phases = 3;         
design.Rm = 0.1;
design.g = 5/1000;
design.Ri = design.Rm + design.g;
design.WmVWp = 0.75;
design.WpVRm = 1.0;
design.RiVRm = design.Ri / design.Rm;
design.RoVRm = 1.2;
design.RaVRo = 1.05;
design.RsoVRm = 0.2;
design.RsiVRso = 0;
design.WcVWp = 1/3;
design.CoilFillFactor = 0.65;

% set the number of turns or the wire diameter or both
%design.Dc = 1/1000;
design.Ntot = 500;
% set the mode of the design
design.mode = 2;
% Ratio of grid resistance to machine resistance
design.RlVRp = 10; 
% ratio of grid inductance to machine inductance
design.LgVLc = 0;
% set the number of Poles in each part to 1 as we will be multiplying up
% the Poles to get the required power specified in optins.targetpower
design.Poles = [1 1];


design = ratios2dimensions_ACTM(design);

RunFEMMSimWithCoils_ACTM(design.WmVWp, design.WpVRm, design.RiVRm, design.RoVRm, design.RsoVRm, design.WcVWp, design.Rm, design.Ntot, design.CoilFillFactor, [0 0 0], design.mode)

% r-component: 0.000221272 Tesla meter^3
% z-component: -3.02541e-008 Tesla meter^3

intBrNoCurrent = 0.000221272;

J = 1e6;

predyForce = intBrNoCurrent * J

% Force via weighted stress tensor
% r-component: 0 N
% z-component: -234.03 N

% B.J
% r-component: 0 N
% z-component: -221.272 N

