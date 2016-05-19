% Test_systemodeforcefcn_linear_mvarm
%
% test with snapper stuff

clear

%% Set up the design

% Common

design.Taup = 83.3/1000;
design.ls = 300/1000;
design.g = 3.5/1000;
nominalg = design.g;
design.ks = 30e3; % 30e3;
design.AngleFromHorizontal = pi/2;


% Translator

design.Taum2 = 34/1000;
design.hm2 = 20/1000;
design.hbi2 = 12.7/1000;

% Armature

design.Taum1 = 22/1000;
design.hm1 = 7/1000;
design.hc = 20/1000;
design.hbi1 = 7.4/1000;
design.Poles(1) = 6;
design.Dco = design.Taup - 5/1000;
design.Dci = (design.Taum1) + (2/1000);
design.CoilFillFactor = 0.6;
design.Dc = 0.75/1000;
design.RlVRp = 15; % 10;
design.extraAMassFact = 1.25;
design.Cd = 0.1; %0; %1.05;
design.DragArea = 0.01;
design.Phases = 1;
design.mu_fA = 0.2;

design = dimensions2ratios_snapper(design);

design.mode = [1, 0, 0, 0];
design.HcMag = [mgoe2hc(35), mgoe2hc(35)];
design.FieldDirection = 1;

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
simoptions.skip = 1;
% the time span of the simulation
simoptions.tspan = [0, 9];   
% Additional absolute tolerances on the components of the solution
simoptions.abstol = [];

simoptions.NoOfMachines = 1;

speed = 0.042;

simoptions.drivetimes = simoptions.tspan(1):((simoptions.tspan(2) - simoptions.tspan(1))/20):simoptions.tspan(2);
simoptions.xT = simoptions.drivetimes .* speed;
simoptions.vT = repmat(speed, size(simoptions.xT));


%% set up the functions and run sim

% objactiam_machine will perform the simulation of each machine, so use a
% dummy simulation function
% simoptions.simfun = 'simfunnocurrent_SNAPPER';
simoptions.simfun = 'dummysimfun';
simoptions.finfun = 'finfun_SNAPPER';
simoptions.odeevfun = 'prescribedmotodeforcefcn_linear_mvgarm';
% simoptions.dpsidxfun = 'polypsidot_ACTIAM'; %@dpsidx_tubular; 
simoptions.resfun = 'prescribedmotresfun_linear_mvgarm';
simoptions.forcefcn = 'forcefcn_linear_mvgarm_pscbmot';
simoptions.ODESim.ForceFcnArgs = {};

[T, Y, results, design] = simulatemachine_linear(design, simoptions, ...
                                                 simoptions.simfun, ...
                                                 simoptions.finfun, ...
                                                 simoptions.odeevfun, ...
                                                 simoptions.resfun);


plotresultsproscribedmot_linear(T, Y, results, design, 1)

%%

[design2, simoptions] = simfunnocurrent_SNAPPER(design, simoptions);
[design2, simoptions] = finfun_SNAPPER(design2, simoptions);

load design_006_wholesys_design_and_simoptions.mat

%%

% plotresultsbuoysys_snapper(T, Y, results, 1)
simoptions.ODESim.InitialConditions = zeros(1, size(Y,2));
simoptions.NoOfMachines = 1;
design.PoleWidth = design.Taup;
[results2, design2] = systemresfun_SNAPPER(T, Y, design2, simoptions);
results2.Fs = results2.Fa(:,1);
results2.Ffea = results2.Fpto + results2.Fa(:,3);
plotresultsbuoysys_snapper(T, Y, results2, 1)

%%

simoptions2 = simoptions;
design2 = design;

% First get the snapper root/trunk directory
snappertrunkdir = fileparts(which('wholesystemsim_Snapper'));

% Set up the sea data file locations, these are for the 2m buoy
simoptions2.HeaveFile = ['Cylinder_2m_dia_d010410/',...
                'heave_coefficients_cyl_2di_1dr_d020610.mat'];
            
simoptions2.SurgeFile = ['Cylinder_2m_dia_d010410/',...
                        'surge_coefficients_cyl_2di_1dr.mat'];                        
                    
simoptions2.HydroCoeffsFile = 'Cylinder_2m_dia_d010410/cyl_d3103v4.1';

simoptions2.ExcitationFile = 'Cylinder_2m_dia_d010410/cyl_d3103v4.2';

simoptions2.BuoyParameters = load(fullfile(getbuoylibdir, ...
                      'Cylinder_2m_dia_d010410',...
                      'buoyparams_d3103v4.mat'));

                  
simoptions2.buoylibdir = getbuoylibdir;

% design2.buoynum = -1;

% objactiam_machine will perform the simulation of each machine, so use a
% dummy simulation function
% if all(isfield(design2, {'p_FEAFy', 'slm_psidot'}))
%     simoptions2.simfun = 'dummysimfun';
%     simoptions2.finfun = 'dummysimfun';
% else
    simoptions2.simfun = 'simfunnocurrent_SNAPPER';
    simoptions2.finfun = 'systemfinfun_SNAPPER';
% end

simoptions2.odeevfun = 'systemodeforcefcn_linear_mvgarm';
% simoptions2.dpsidxfun = 'polypsidot_ACTIAM'; %@dpsidx_tubular; 
simoptions2.resfun = 'systemresfun_SNAPPER';
simoptions2.preallocresfcn = 'prallocresfcn_SNAPPER';
simoptions2.forcefcn = 'forcefcn_snapper';
simoptions2.ODESim.ForceFcnArgs = {};

if isfield(design2, 'Ntot')
    design2 = rmfield(design2, 'Ntot');
end

[T2, Y2, results2, design2, simoptions2] = simulatemachine_linear(design2, ...
                                                 simoptions2, ...
                                                 simoptions2.simfun, ...
                                                 simoptions2.finfun, ...
                                                 simoptions2.odeevfun, ...
                                                 simoptions2.resfun);

results2.Fs = results2.Fa(:,1);
results2.Ffea = results2.Fpto + results2.Fa(:,3);
plotresultsbuoysys_snapper(T2, Y2, results2, 1)



