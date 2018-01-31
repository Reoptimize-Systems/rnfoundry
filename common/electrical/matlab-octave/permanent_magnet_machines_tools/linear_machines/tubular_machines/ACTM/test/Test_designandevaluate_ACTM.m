%% Test_designAndEvaluate_ACTM

clear 
% First we need some plausible machine variables for testing, we will use
% the AWS variables for this

design.Phases = 3;         % Number of Phases in machine
design.Rm = 0.1;
design.g = 5/1000;
design.Ri = design.Rm + design.g;
design.WmVWp = 0.75;
design.WpVRm = 0.5;
design.RiVRm = design.Ri / design.Rm;
design.RoVRm = 1.2;
design.RaVRo = 1.025;
design.RsoVRm = 0.1;
design.RsiVRso = 0;
design.WcVWp = 1/3;
design.CoilFillFactor = 0.65;
%design.Dc = 1/1000;  % 1 mm diameter wire 
design.Ntot = 500;
design.mode = 2; 
design.LgVLc = 0;
design.Poles = [1 1];
% Acon = pi * (dc/2)^2;
% Acu = (Rm * WpVRm * WcVWp) * ((RoVRm * Rm) - Ri) * kfill;
% AconVAcu = Acon / Acu;

% E = [200e9 151e9];  % Young's modulus of elasticity for Structural steel
% and laminated steel respectively

design = ratios2dimensions_ACTM(design);

options.E = [207e9 100e5 207e9];
simoptions.Evaluation.targetPower = 10e3; % 10kW machine
simoptions.Evaluation.mlength = 4; % Overlap between stator and translator, i.e. stator is mleng metres longer than the translator
% options.pointsPerPole = 40;
options.coilYieldStrength = 70e6;

% Test with linear motion
speed = 1;
simoptions.ODESim.InitialConditions = zeros(1, design.Phases);
simoptions.ODESim.ResultsTSkip = 1;
simoptions.ODESim.TimeSpan = [0, 5];
simoptions.drivetimes = 0:simoptions.ODESim.TimeSpan(2);
simoptions.vT = repmat(speed, size(simoptions.drivetimes));
simoptions.xT = simoptions.vT .* simoptions.drivetimes;
simoptions.Lmode = 1;

% set up the functions
simoptions.ODESim.PreProcFcn = 'simfun_ACTM';
simoptions.ODESim.PostPreProcFcn = 'prescribedmotfinfun_ACTM';
simoptions.ODESim.EvalFcn = 'prescribedmotodeforcefcn_linear';
simoptions.ODESim.ForceFcn = 'forcefcn_linear_pscbmot';
simoptions.ODESim.PostSimFcn = 'prescribedmotresfun_linear';

%%

design.RlVRp = 10; % Ratio of grid resistance to machine resistance

[score, design, simoptions, T, Y, results] = designandevaluate_ACTM(design, simoptions);

%%
% if ~isfemmopen
%     openfemm;
%     %main_minimize;
% end
% 
% RunFEMMSimWithCoils_ACTM_Dc_Value = RunFEMMSimWithCoils_ACTM(design.WmVWp, design1.WpVRm, design1.RiVRm, design1.RoVRm, design1.RsoVRm, design1.WcVWp, design1.Rm, design1.Ntot, design1.CoilFillFactor, [0 0 0], [0 1]);
% 




