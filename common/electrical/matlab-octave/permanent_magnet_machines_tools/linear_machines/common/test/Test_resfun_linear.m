% Test_resfun_linear
                    

%% ACTM

clear design

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
% design.FieldDirection = 1;
% design.PowerPoles = Poles(1);

design = ratios2dimensions_ACTM(design);

simoptions.simfun = 'simfun_ACTM';
mname = 'ACTM';

%% Set up Common Parameters

design.RlVRp = 10;

simoptions.Lmode = 0;
simoptions.NoOfMachines = 1;

%% Test with linear motion

speed = 1;
simoptions.IC = zeros(1, design.Phases);
simoptions.skip = 1;
simoptions.tspan = [0, 5];
simoptions.drivetimes = 0:simoptions.tspan(2)/2:simoptions.tspan(2);
simoptions.vT = repmat(speed, size(simoptions.drivetimes));
simoptions.xT = simoptions.vT .* simoptions.drivetimes;
simoptions.tether_length = 0;
simoptions.NoOfMachines = 1;

simoptions.odeevfun = 'simplelinearmachineode_proscribedmotion'; 
simoptions.resfun = 'resfun_linear';
simoptions.finfun = ['finfun_', mname];

[T, Y, results, design] = simulatemachine_linear(design, simoptions, simoptions.simfun, ...
                                                  simoptions.finfun, simoptions.odeevfun, simoptions.resfun); 

% plotresultsproscribedmot_linear(T, Y, results, 1);                                              

%% Test with sinusoidal motion

simoptions.IC = zeros(1, design.Phases);
simoptions.skip = 1;
simoptions.xTperiod = 3;
simoptions.xTamplitude = 1;
simoptions.tspan = [0, 9];
simoptions.drivetimes = 0:simoptions.tspan(2)/100:simoptions.tspan(2);
simoptions.xT = simoptions.xTamplitude * sin(2 .* pi .* (1/simoptions.xTperiod) .* simoptions.drivetimes - pi/2);
simoptions.vT = 2 .* pi .* (1/simoptions.xTperiod) .* simoptions.xTamplitude * cos(2 .* pi .* (1/simoptions.xTperiod) .* simoptions.drivetimes - (pi/2));
simoptions.Lmode = 0;
simoptions.tether_length = 0;
simoptions.NoOfMachines = 1;

simoptions.odeevfun = 'simplelinearmachineode_proscribedmotion'; 
simoptions.resfun = 'resfun_linear';
simoptions.finfun = ['finfun_', mname];

[T, Y, results, design] = simulatemachine_linear(design, simoptions, simoptions.simfun, ...
                                                  simoptions.finfun, simoptions.odeevfun, simoptions.resfun); 
                                              
plotresultsproscribedmot_linear(T, Y, results, 1);

%% Test with buoy in sinusoidal sea

% Set up the buoy and sea data files, these are for the 2m buoy
% snappertrunkdir = fileparts(which('wholesystemsim_Snapper'));
simoptions.HeaveFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410', 'heave_coefficients_cyl_2di_1dr_d020610.mat');
simoptions.SurgeFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410', 'surge_coefficients_cyl_2di_1dr.mat');
simoptions.HydroCoeffsFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410','cyl_d3103v4.1');
simoptions.ExcitationFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410','cyl_d3103v4.2');
simoptions.BuoyParameters = load(fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410', 'buoyparams_d3103v4.mat'));
simoptions.buoy = [];

simoptions.tspan = [0, 60];
% params.amp = 1;
params.sigma = 2 * pi * 0.35;
params.phase = pi/2;

simoptions.SeaParameters = defaultseaparamaters(params);

simoptions.tether_length = 4;

simoptions.odeevfun = 'systemode_linear'; 
simoptions.finfun = ['systemfinfun_', mname];
simoptions.resfun = 'systemresfun_linear';    

[T, Y, results, design] = simulatemachine_linear(design, simoptions, simoptions.simfun, ...
                                                  simoptions.finfun, simoptions.odeevfun, simoptions.resfun); 

% plotresultsbuoysys_linear(T, Y, results, design, 1)

