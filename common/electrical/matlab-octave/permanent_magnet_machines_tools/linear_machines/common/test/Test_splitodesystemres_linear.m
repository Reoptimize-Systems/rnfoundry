%% Test_splitodesystemres_linear

clear 
% First we need some plausible machine variables for testing, we will use
% the AWS variables for this

design.phases = 3;         % Number of phases in machine
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
design.fillfactor = 0.65;
%design.Dc = 1/1000;  % 1 mm diameter wire 
design.Ntot = 500;
design.mode = 2; 
design.LgVLc = 0;
design.poles = [10 30];
design.RgVRc = 10;

design = ratios2dimensions_ACTM(design);

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
% determined by a penalty. The number of poles in the design will be
% modified if exceeded
simoptions.maxAllowedTLength = 5;
simoptions.maxAllowedxT = inf;
% determines method used to calculate inductance
simoptions.Lmode = 1;
% the initial values of xA, vA and the initial currents in the coils at t=0
simoptions.IC = [0, 0, zeros(1, design.phases)];
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
simoptions.buoynum = -1;

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

%% set up the functions and run sim

% objactiam_machine will perform the simulation of each machine, so use a
% dummy simulation function
% simoptions.simfun = 'dummysimfun';
simoptions.simfun = 'systemsimfun_ACTM';
simoptions.finfun = 'systemfinfun_ACTM';
simoptions.odeevfun = 'systemode_linear';
% simoptions.dpsidxfun = 'polypsidot_ACTIAM'; %@dpsidx_tubular; 
% simoptions.resfun = 'systemresfun_linear_mvgarm';
simoptions.resfun = 'splitsystemresfun_linear';
simoptions.spfcn = 'splitodesystemres_linear';
simoptions.splitode = 10;
simoptions.preallocresfcn = 'springandsnapprallocresfcn_linear_mvgarm';
% simoptions.forcefcn = 'springforcefcn_linear_mvgarm';
% simoptions.forcefcnargs = {};
% simoptions.forcefcn = 'springandsnapforcefcn_linear_mvgarm';
% simoptions.forcefcnargs = {};

[T, Y, results, design] = simulatemachine_linear(design, simoptions, ...
                                                 simoptions.simfun, ...
                                                 simoptions.finfun, ...
                                                 simoptions.odeevfun, ...
                                                 simoptions.resfun);

                                             
%% Compare to straight sim

% objactiam_machine will perform the simulation of each machine, so use a
% dummy simulation function
% simoptions.simfun = 'dummysimfun';
simoptions.simfun = 'systemsimfun_ACTM';
simoptions.finfun = 'systemfinfun_ACTM';
simoptions.odeevfun = 'systemode_linear';
% simoptions.dpsidxfun = 'polypsidot_ACTIAM'; %@dpsidx_tubular; 
% simoptions.resfun = 'systemresfun_linear_mvgarm';
simoptions.resfun = 'systemresfun_linear';

if isfield(simoptions, 'spfcn')
    simoptions = rmfield(simoptions, 'spfcn');
end

if isfield(simoptions, 'splitode')
    simoptions = rmfield(simoptions, 'splitode');
end

simoptions.preallocresfcn = 'prallocresfcn_linear_system';

% simoptions.forcefcn = 'springforcefcn_linear_mvgarm';
% simoptions.forcefcnargs = {};
% simoptions.forcefcn = 'springandsnapforcefcn_linear_mvgarm';
% simoptions.forcefcnargs = {};

[T, Y, results2, design2] = simulatemachine_linear(design, simoptions, ...
                                                 simoptions.simfun, ...
                                                 simoptions.finfun, ...
                                                 simoptions.odeevfun, ...
                                                 simoptions.resfun);                                             
