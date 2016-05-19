% Test_systemodeforcefcn_linear_mvarm_ACTM
%
% test with snapper stuff

clear

%% Set up the design

load actm_design_and_simoptions.mat

clear simoptions

design.Poles = [15, 4];

design.PowerPoles = design.Poles(2);

design.ks = 30000;

design.massA = design.Phases * design.Poles(2) * design.wlength  * (pi * (design.Dc/2).^2) * 8600;

design.weightA = 9.81 * design.massA;

design.massT = fieldpoleweight_TM(design.WmVWp, design.WpVRm, ...
                              design.RsiVRso, design.RsoVRm, ...
                              design.Rm, 7500, 7500, 7500) * design.Poles(1) / 9.81;

% Determine the amount the spring is compressed by the weight of
% the armature. It will start from this compressed position
design.xSC = design.weightA / design.ks;

simoptions.maxAllowedxA = 1e6;
design.maxAllowedxA = simoptions.maxAllowedxA;

design.sides = 1;
design.Phases = 3;

designACTM = design;


designACTM.Dc = 1/1000;
[designACTM.Ntot, designACTM.Dc] = CoilTurns(designACTM.Hc * designACTM.Wc, design.CoilFillFactor, designACTM.Dc);
designACTM.CoilResistance = designACTM.wlength * 1.68e-8  / (pi*(designACTM.Dc/2)^2);
designACTM.LoadResistance = designACTM.CoilResistance * 10;
designACTM.R = designACTM.CoilResistance + designACTM.LoadResistance;
designACTM.L = designACTM.L .* designACTM.Ntot / design.Ntot;


load design_006.mat;

% Make snapping force of same size as original machine
designACTM.slm_Fsnap = slmengine(mdata.indepvar(45:55,1), mdata.FEAFy(45:55,1) .* design.ls .* 6,...
        'plot','off',...
        'verbosity',0,...
        'knots',6,...
        'leftvalue', 0,...
        'rightvalue', 0,...
        'InteriorKnots', 'free');

%% Set up the simulation options

% set up the penalties, empty penalties will be ignored

% the maximum allowed rms current density in a coil
simoptions.maxAllowedJrms = 6e6;
% the maximum allowed peak current density in the coil
simoptions.maxAllowedJpeak = 10e6;
% the minimum allowed rms voltage produced in a coil
simoptions.minAllowedRMSEMF = 200;
% the maximum allowed voltage produced by a coil
simoptions.maxAllowedEMFpeak = [];

% set the other simulation parameters 

% The maximum allowed translator length, this is a hard limit, not
% determined by a penalty. The number of Poles in the design will be
% modified if exceeded
simoptions.maxAllowedTLength = 5;
% determines method used to calculate inductance
simoptions.Lmode = 1;
% the initial values of xA, vA and the initial currents in the coils at t=0
simoptions.ODESim.InitialConditions = [0, 0, zeros(1, designACTM.Phases)];
% the number of calculations to skip when producing output after the ode
% solver finishes
simoptions.skip = 1;
% the time span of the simulation
simoptions.tspan = [0, 60];   
% Additional absolute tolerances on the components of the solution
simoptions.abstol = [];
% The number of machines attached to the buoy
simoptions.NoOfMachines = 1;

%% Set up the sea and buoy - prototype buoy vals

% First get the snapper root/trunk directory
snappertrunkdir = fileparts(which('wholesystemsim_Snapper'));
% Set up the sea data file locations, these are for the 2m buoy
simoptions.HeaveFile = fullfile('Cylinder_2m_dia_d010410',...
                'heave_coefficients_cyl_2di_1dr_d020610.mat');
            
simoptions.SurgeFile = fullfile('Cylinder_2m_dia_d010410',...
                        'surge_coefficients_cyl_2di_1dr.mat');     
                    
simoptions.HydroCoeffsFile = fullfile('Cylinder_2m_dia_d010410', 'cyl_d3103v4.1');

simoptions.ExcitationFile = fullfile('Cylinder_2m_dia_d010410','cyl_d3103v4.2');

simoptions.BuoyParameters = load(fullfile(getbuoylibdir, ...
                      'Cylinder_2m_dia_d010410',...
                      'buoyparams_d3103v4.mat'));

simoptions.buoylibdir = getbuoylibdir;
simoptions.buoy = [];

% Use a random sea with a peak frequency of 0.35 Hz
params.peak_freq = 0.35;
params.phase = pi./ 2;
% Call defaultseaparamaters to set up the necessary sea data
simoptions.SeaParameters = defaultseaparamaters(params);

% % Use a random sea with a peak frequency of 0.35 Hz
% params.peak_freq = 0.35;
% params.phase = pi./2;
% % Call defaultseaparamaters to set up the necessary sea data
% simoptions.SeaParameters = defaultseaparamaters(params);

% set the initial tether length between the buoy and the hawser
simoptions.tether_length = 3;

%% Set up the sea and buoy - larger buoy

% The peak frequency of a PM Specturm
SeaParameters.peak_freq = 1/9;
% SeaParameters.peak_freq = 0.35;
SeaParameters.phase = pi/2;
% % The range of frequencies to be included in the spectrum
SeaParameters.sigma_range = [2*pi*0.055, 0.1*(0.6)^0.5 + 2.24];
% % The number of frequencies to be produced
% SeaParameters.freqcount = 50;
% % The depth of the water
SeaParameters.water_depth = 50;

simoptions.SeaParameters = defaultseaparamaters(SeaParameters);

% 
simoptions.tether_length = 6;
% 
% simoptions.FarmSize = 10e6;
% 
% use buoy number 37, 4m diameter, 2m draft
simoptions.buoynum = 37;
designACTM.buoynum = simoptions.buoynum;

%% set up the functions and run sim

% objactiam_machine will perform the simulation of each machine, so use a
% dummy simulation function
simoptions.simfun = 'dummysimfun';
simoptions.finfun = 'systemfinfun_ACTM';
simoptions.odeevfun = 'systemodeforcefcn_linear_mvgarm';
% simoptions.dpsidxfun = 'polypsidot_ACTIAM'; %@dpsidx_tubular; 
% simoptions.resfun = 'systemresfun_linear_mvgarm';
simoptions.resfun = 'systemresfun_linear_mvgarm';
simoptions.preallocresfcn = 'springandsnapprallocresfcn_linear_mvgarm';
% simoptions.forcefcn = 'springforcefcn_linear_mvgarm';
% simoptions.ODESim.ForceFcnArgs = {};
simoptions.forcefcn = 'springandsnapforcefcn_linear_mvgarm';
simoptions.ODESim.ForceFcnArgs = {};

[T, Y, results, designACTM] = simulatemachine_linear(designACTM, simoptions, ...
                                                 simoptions.simfun, ...
                                                 simoptions.finfun, ...
                                                 simoptions.odeevfun, ...
                                                 simoptions.resfun);

results.Ffea = results.Fpto + results.Fa(:,2);
results.Fs = results.Fa(:,1);
plotresultsbuoysys_snapper(T, Y, results, 1)
% figure; plot(T, results.Fpto, T, results.Fa(:,2));
figure; plot(T, results.Fpto);
addaxis(T, results.vT - Y(:,6));


