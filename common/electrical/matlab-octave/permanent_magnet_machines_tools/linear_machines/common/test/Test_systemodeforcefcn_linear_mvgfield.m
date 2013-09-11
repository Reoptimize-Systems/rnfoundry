% Test_systemodeforcefcn_linear_mvgfield
%
% test with snapper stuff

clear; close all;

%% Set up the design

design.bpVTaup = 0.85; 
design.lmVbp = 0.2; 
design.dgVlm = 1.4; 
design.lsVTaup = 3;
design.dbiVlm = 0.3;
design.WcVTaup = 1/3;
design.hcVgap = 0.95;
design.Taup = 0.1;
design.fillfactor = 0.55;
design.J = 0;
design.klineardrag = 0;
design.mu_fA = 0;
design.mu_fF = 0;
% design.massF = 100;
% design.weightF = 9.81 * design.massF;
design.mu_fT = 0;
design.massT = 0;
design.xSC = 0;
design.maxAllowedxA = inf;

design.Cd = 0;
design.DragArea = 0;


design = ratios2dimensions_ACPMSM(design);

design.AngleFromHorizontal = pi/2;
design.phases = 3;
design.sides = 2; %hmmm
design.poles = [12 15];
design.Ntot = 500;
design.RgVRc = 10;
design.LgVLc = 0;
design.OuterWebs = 3;
design.GuideRailIMethod = '1.3';
design.GuideRailIVars = [0.1, 0.1, 0.095, 0.095];
% 
% % Get the dimensionless ratios from the parameters
% design = dimensions2ratios_ACPMSM(design);

% % set design mode
% design.mode = [1, 1, 0, 1];

design.OuterWebs = 3;
design.GuideRailIMethod = '1.3';
design.GuideRailIVars = [0.1, 0.1, 0.095, 0.095];
design.InnerStructureBeamVars = [];

% Calculate the extra length needed for the bearings
bearingWidth = 0.1; 
options.alphab = (design.ls + 2*bearingWidth) / design.ls;

% Common

design.ks = 4 * 60e3;
% design.FieldDirection = -1;

% design.HcMag = mgoe2hc(35);

design.MagCouple.hm = 0.025;
design.MagCouple.hms = 0.01;
design.MagCouple.g = 5/1000;
design.MagCouple.ht = 0.1;
design.MagCouple.htbi = design.MagCouple.ht/5;

% MagCouple.Wt = 1*MagCouple.ht;
design.MagCouple.Wt = 0.05;
% MagCouple.Wt = 0.015;
% MagCouple.Wr = 5*MagCouple.Wt;
% MagCouple.Wr = 0.25;
% MagCouple.Wr = 0.05;
design.MagCouple.Wr = 0.1;

design.MagCouple.Wms = 0.95*design.MagCouple.Wr;

design.MagCouple.Wm = 1.2 * design.MagCouple.Wt;

design.MagCouple.ls = 0.2;

design.MagCouple.FieldSteelDensity = 7500;

design.MagCouple.FieldMagDensity = 7500;

design.MagCouple.N = 8;

%%

simoptions = buoysimoptions;

% speed = 1;
% simoptions.IC = 0;
% simoptions.skip = 1;
% simoptions.tspan = [0, 10];
% simoptions.drivetimes = 0:simoptions.tspan(2);
% simoptions.vT = repmat(speed, size(simoptions.drivetimes));
% simoptions.xT = simoptions.vT .* simoptions.drivetimes;
% simoptions.Lmode = 1;

% use buoy number 37, 4m diameter, 2m draft
simoptions.buoynum = 37;

simoptions = buoysetup(simoptions.buoynum, [], simoptions);

simoptions.SeaParameters = seasetup('MinFreq', simoptions.BuoyParameters.minfreq, ...
                                    'MaxFreq', simoptions.BuoyParameters.maxfreq, ...
                                    'PMPeakFreq', 0.111, ...
                                    'NoOfFrequencies', 50, ...
                                    'WaterDepth', simoptions.BuoyParameters.water_depth);


%simoptions.SeaParameters.sigma = 2 * pi * 0.111;
%simoptions.SeaParameters.phase = pi / 2;
% simoptions.SeaParameters.peak_freq = 0.111;
% simoptions.SeaParameters.water_depth = 40;
% simoptions.SeaParameters = defaultseaparamaters(simoptions.SeaParameters);

simoptions.tether_length = 5;

simoptions.NoOfMachines = 1;

simoptions.tspan = [0, 120];
simoptions.Lmode = 1;
simoptions.IC = [0, 0, 0];

simoptions.simfun = 'simfun_MAGCOUPLE';
simoptions.finfun = 'finfun_MAGCOUPLE';
% simoptions.odeevfun = @simplelinearmachineode_proscribedmotion;
simoptions.odeevfun = 'systemodeforcefcn_linear_mvgfield';
simoptions.resfun = 'systemresfun_linear_mvgfield';
% simoptions.dpsidxfun = @polydpsidx_ACPMSM
simoptions.forcefcn = 'forcefcn_linear_mvgfield_system';
simoptions.forcefcnargs = {};

simfunargs = {@simfun_ACPMSM};

finfunargs = {@systemfinfun_ACPMSM};

simoptions.rho = 0;

% design.FieldDirection = -1;

%%

[T, Y, results, design] = simulatemachine_linear(design, simoptions, ...
                                                 simoptions.simfun, ...
                                                 simoptions.finfun, ...
                                                 simoptions.odeevfun, ...
                                                 simoptions.resfun, ...
                                                 'simfunargs', simfunargs, ...
                                                 'finfunargs', finfunargs);
                                             
                                             
plotresultsbuoysys_linear_mvgfield(T, Y, results, design, 1)

design.GridMeanPower

machinescore_ACPMSM()


% %% no FEA
% 
% simoptions.simfun = 'dummysimfun';
% simfunargs = {@dummysimfun};
% 
% 
% design.MagCouple.N = 8;
% design.Ntot = 1000;
% 
% 
% [T, Y, results, design] = simulatemachine_linear(design, simoptions, ...
%                                                  simoptions.simfun, ...
%                                                  simoptions.finfun, ...
%                                                  simoptions.odeevfun, ...
%                                                  simoptions.resfun, ...
%                                                  'simfunargs', simfunargs, ...
%                                                  'finfunargs', finfunargs);
%                                              
%                                              
% plotresultsbuoysys_linear_mvgfield(T, Y, results, design, 1)
% 
% design.GridMeanPower


%% Compare to no magcouple


simoptions.simfun = 'dummysimfun';
% simoptions.simfun = 'simfun_ACPMSM';
simoptions.finfun = 'systemfinfun_ACPMSM';
simoptions.odeevfun = 'systemode_linear';
simoptions.resfun = 'systemresfun_linear';

simoptions.rho = 0;

[T, Y, results, design] = simulatemachine_linear(design, simoptions, ...
                                                 simoptions.simfun, ...
                                                 simoptions.finfun, ...
                                                 simoptions.odeevfun, ...
                                                 simoptions.resfun);
                                             
                                             
plotresultsbuoysys_linear(T, Y, results, design, 1)


design.GridMeanPower


