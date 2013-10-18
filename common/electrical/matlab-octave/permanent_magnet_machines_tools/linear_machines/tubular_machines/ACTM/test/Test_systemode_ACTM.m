%% Test_systemode_ACTM

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
design.Poles = [10 30];
% Acon = pi * (dc/2)^2;
% Acu = (Rm * WpVRm * WcVWp) * ((RoVRm * Rm) - Ri) * kfill;
% AconVAcu = Acon / Acu;

% E = [200e9 151e9];  % Young's modulus of elasticity for Structural steel
% and laminated steel respectively

options.E = [207e9 100e5 207e9];
options.targetPower = 10e3; % 10kW machine
options.mlength = 4; % Overlap between stator and translator, i.e. stator is mleng metres longer than the translator
% options.pointsPerPole = 40;
options.coilYieldStrength = 70e6;

% Chose the method used to evaluate inductance
simoptions.Lmode = 1;

% Test with buoy system
% Set up the buoy and sea data files, these are for the 2m buoy
% snappertrunkdir = fileparts(which('wholesystemsim_Snapper'));
simoptions.HeaveFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410', 'heave_coefficients_cyl_2di_1dr_d020610.mat');
simoptions.SurgeFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410', 'surge_coefficients_cyl_2di_1dr.mat');
simoptions.HydroCoeffsFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410','cyl_d3103v4.1');
simoptions.ExcitationFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410','cyl_d3103v4.2');
simoptions.BuoyParameters = load(fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410', 'buoyparams_d3103v4.mat'));

params.amp = 1;
params.peak_freq = 0.35; % centred at resonant frequency
params.phase = pi/2;

simoptions.SeaParameters = defaultseaparamaters(params);

simoptions.tether_length = 4;

simoptions.tspan = [0, 30];

% set up the functions
simoptions.simfun = 'simfun_ACTM';
simoptions.finfun = 'systemfinfun_ACTM';
simoptions.odeevfun = 'systemode_linear'; 
simoptions.resfun = 'systemresfun_linear';

%%

design.RlVRp = 10; % Ratio of grid resistance to machine resistance

[score, design, simoptions, T, Y, results] = designandevaluate_ACTM(design, simoptions, options);

%%

% plotresultsbuoysys_linear(T, Y, results, 1)

%%

% fps = 30;
% tscale = 1;
% frames = tscale*round((max(T)-min(T))*fps);
% animatesytem_TM(design, simoptions, T, Y, results, frames, fps, 'supergen_ACTM_mono_waves.avi')

%%
% if ~isfemmopen
%     openfemm;
%     %main_minimize;
% end
% 
% RunFEMMSimWithCoils_ACTM_Dc_Value = RunFEMMSimWithCoils_ACTM(design.WmVWp, design1.WpVRm, design1.RiVRm, design1.RoVRm, design1.RsoVRm, design1.WcVWp, design1.Rm, design1.Ntot, design1.CoilFillFactor, [0 0 0], [0 1]);
% 




