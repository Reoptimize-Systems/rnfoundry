% Test_systemodeforcefcn_linear_mvgarm_snapper_magcouple
%
%
% test with snapper stuff

close all; clear

%% Set up the design

% Common

design.Taup = 83.3/1000;
design.ls = 300/1000;
design.g = 3/1000;
nominalg = design.g;
design.ks = 30e3;


% Translator

design.Taum2 = 35/1000;
design.hm2 = 20/1000;
design.hbi2 = 12.7/1000;

% Armature

design.Taum1 = 22/1000;
design.hm1 = 5/1000;
design.hc = 20/1000;
design.hbi1 = 7.4/1000;
design.Poles(1) = 6;
design.Dco = design.Taup - 5/1000;
design.Dci = (design.Taum1) + (2/1000);
design.CoilFillFactor = 0.6;
design.Dc = 0.75/1000;
design.RlVRp = 10;
design.extraAMassFact = 0.25;
design.Cd = 1.05;
design.DragArea = 0.01;
design.Phases = 1;

design = dimensions2ratios_snapper(design);

design.mode = [1, 0, 0, 0];
design.HcMag = mgoe2hc(35);




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
simoptions.ODESim.InitialConditions = [0, 0, 0];
% the number of calculations to skip when producing output after the ode
% solver finishes
simoptions.ODESim.ResultsTSkip = 1;
% the time span of the simulation
simoptions.ODESim.TimeSpan = [0, 120];   
% Additional absolute tolerances on the components of the solution
simoptions.abstol = [];

simoptions.buoynum = -1;
design.buoynum = simoptions.buoynum;

% First get the snapper root/trunk directory
% snappertrunkdir = fileparts(which('wholesystemsim_Snapper'));
% Set up the sea data file locations, these are for the 2m buoy
simoptions.HeaveFile = fullfile('Cylinder_2m_dia_d010410',...
                'heave_coefficients_cyl_2di_1dr_d020610.mat');
            
simoptions.SurgeFile = fullfile('Cylinder_2m_dia_d010410',...
                        'surge_coefficients_cyl_2di_1dr.mat');     
                    
simoptions.HydroCoeffsFile = fullfile('Cylinder_2m_dia_d010410', 'cyl_d3103v4.1');

simoptions.ExcitationFile = fullfile('Cylinder_2m_dia_d010410','cyl_d3103v4.2');

simoptions.BuoySim.BuoyParameters = load(fullfile(getbuoylibdir, ...
                      'Cylinder_2m_dia_d010410',...
                      'buoyparams_d3103v4.mat'));

simoptions.buoylibdir = getbuoylibdir;


% Use a sinusoidal sea with a frequency of 0.35 Hz
params.sigma = 2 * pi * 0.35;
params.phase = pi/2;
% Call defaultseaparamaters to set up the necessary sea data
simoptions.BuoySim.SeaParameters = defaultseaparamaters(params);

% % Use a random sea with a peak frequency of 0.35 Hz
% params.peak_freq = 0.35;
% % Call defaultseaparamaters to set up the necessary sea data
% simoptions.BuoySim.SeaParameters = defaultseaparamaters(params, simoptions.BuoySim.BuoyParameters);

% set the initial tether length between the buoy and the hawser
simoptions.BuoySim.tether_length = 3;

% % The peak frequency of a PM Specturm
% SeaParameters.peak_freq = 1/9;
% % The range of frequencies to be included in the spectrum
% SeaParameters.sigma_range = [2*pi*0.055, 0.1*(0.6)^0.5 + 2.24];
% % The number of frequencies to be produced
% SeaParameters.freqcount = 50;
% % The depth of the water
% SeaParameters.water_depth = 50;
% 
% simoptions.BuoySim.SeaParameters = defaultseaparamaters(SeaParameters);
% 
% simoptions.BuoySim.tether_length = 10;
% 
% simoptions.FarmSize = 10e6;
% 
% use buoy number 37, 4m diameter, 2m draft
% simoptions.buoynum = 37;
% design.buoynum = simoptions.buoynum;

simoptions.NoOfMachines = 1;

simoptions.mx_initial_conditions = [0, 0, 0];



%% set up the functions and run sim

[design, simoptions] = simfunnocurrent_SNAPPER(design, simoptions);

[design, simoptions] = systemfinfun_SNAPPER(design, simoptions);


% set up the magnetic coupling with appropriate values to replicate the
% original snapper configuration
design.MagCouple.hm = [];
design.MagCouple.hms = [];
design.MagCouple.g = [];
design.MagCouple.ht = [];
design.MagCouple.htbi = [];
design.MagCouple.Wt = [];
design.MagCouple.Wms = [];
design.MagCouple.Wm = [];
design.MagCouple.ls = [];

design.MagCouple.Wr = design.Taup;

design.MagCouple.slm_FEAFy = design.slm_FEAFy;

design.MagCouple.N = design.PowerPoles * design.sides;


% objactiam_machine will perform the simulation of each machine, so use a
% dummy simulation function
simoptions.ODESim.PreProcFcn = 'dummysimfun';
simoptions.ODESim.PostPreProcFcn = 'dummysimfun';
simoptions.ODESim.EvalFcn = 'systemodeforcefcn_linear_mvgarm';
% simoptions.dpsidxfun = 'polypsidot_ACTIAM'; %@dpsidx_tubular; 
% simoptions.ODESim.PostSimFcn = 'systemresfun_SNAPPER';
simoptions.ODESim.PostSimFcn = 'systemresfun_linear_mvgarm';
simoptions.preallocresfcn = 'prallocresfcn';
simoptions.ODESim.ForceFcn = 'forcefcn_linear_mvgarm_system';
simoptions.ODESim.ForceFcnArgs = {};

simoptions.rho = 0;

[T, Y, results, design] = simulatemachine_linear(design, simoptions, ...
                                                 simoptions.ODESim.PreProcFcn, ...
                                                 simoptions.ODESim.PostPreProcFcn, ...
                                                 simoptions.ODESim.EvalFcn, ...
                                                 simoptions.ODESim.PostSimFcn);
                                             
plotresultsbuoysys_linear_mvgarm(T, Y, results, design, 1)   

%% compare to new

design.MagCouple = rmfield(design.MagCouple, 'slm_FEAFy');
design = rmfield(design, 'slm_FEAFy');

doplot = true;
dosim_snapper2(design, simoptions, doplot);
