% Test_simulatemachine_0

%% PMSM

clear design

design.Phases = 3;
design.Wp = 0.12;
design.Wm = 0.8*design.Wp;
design.hm = 0.015;
design.kw = 0.84;
%Ns = 470;
%Ns = 6; % ???
design.Hc = 979000;
%ht = 0.1;
design.CoilFillFactor = 0.585;
design.g = 0.003; 
design.ls = 0.2; 
design.Dc = 0.005; 
design.E = [200e9 151e9];
design.Wc=design.Wp/design.Phases;
design.Ws=design.Wc/2; 
design.Wt=design.Wc/2;
design.ht=5*design.Wt;
design.hbf = design.hm;
design.hba = design.hbf;


design.Poles = [15 5];
design.RlVRp = 10;
design.LgVLc = 0;

% Get the dimensionless ratios from the parameters
design = dimensions2ratios_PMSM(design);

% set design mode
design.mode = [0, 1, 0, 1];

% set up the functions
simoptions.ODESim.PreProcFcn = 'simfunnocurrent_PMSM';
mname = 'PMSM';
                                      
%% ACPMSM

clear design 

design.bpVTaup = 0.85; 
design.lmVbp = 0.2; 
design.dgVlm = 2.0; 
design.lsVTaup = 3;
design.dbiVlm = 1;
design.WcVTaup = 1/3;
design.hcVgap = 0.95;
design.Taup = 0.2;
design.Ntot = 1000;
design.CoilFillFactor = 0.55;
design.J = 0;

design = ratios2dimensions_ACPMSM(design);

design.Phases = 3;
design.Poles = [10 30];
design.Ntot = 1000;
design.RlVRp = 10;
design.LgVLc = 0;

% set up the functions
simoptions.ODESim.PreProcFcn = 'simfun_ACPMSM';
mname = 'ACPMSM';

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

simoptions.ODESim.PreProcFcn = 'simfun_ACTM';
mname = 'ACTM';

%% Set up Common Parameters

design.RlVRp = 10;

simoptions.Lmode = 0;
simoptions.NoOfMachines = 1;
simoptions.maxAllowedxT = 0.5;

%% Test with linear motion

speed = 1;
simoptions.ODESim.InitialConditions = zeros(1, design.Phases);
simoptions.ODESim.ResultsTSkip = 1;
simoptions.ODESim.TimeSpan = [0, 5];
simoptions.drivetimes = 0:simoptions.ODESim.TimeSpan(2)/2:simoptions.ODESim.TimeSpan(2);
simoptions.vT = repmat(speed, size(simoptions.drivetimes));
simoptions.xT = simoptions.vT .* simoptions.drivetimes;
simoptions.BuoySim.tether_length = 0;
simoptions.NoOfMachines = 1;

simoptions.ODESim.EvalFcn = 'simplelinearmachineode_proscribedmotion'; 
simoptions.ODESim.PostSimFcn = 'resfun_linear_pscbmot';
simoptions.ODESim.PostPreProcFcn = ['finfun_', mname];

[T, Y, results, design] = simulatemachine_linear(design, simoptions, simoptions.ODESim.PreProcFcn, ...
                                                  simoptions.ODESim.PostPreProcFcn, simoptions.ODESim.EvalFcn, simoptions.ODESim.PostSimFcn); 

plotresultsproscribedmot_linear(T, Y, results, 1);                                              

%% Test with sinusoidal motion

simoptions.ODESim.InitialConditions = zeros(1, design.Phases);
simoptions.ODESim.ResultsTSkip = 1;
simoptions.xTperiod = 3;
simoptions.xTamplitude = 1;
simoptions.ODESim.TimeSpan = [0, 9];
simoptions.drivetimes = 0:simoptions.ODESim.TimeSpan(2)/100:simoptions.ODESim.TimeSpan(2);
simoptions.xT = simoptions.xTamplitude * sin(2 .* pi .* (1/simoptions.xTperiod) .* simoptions.drivetimes - pi/2);
simoptions.vT = 2 .* pi .* (1/simoptions.xTperiod) .* simoptions.xTamplitude * cos(2 .* pi .* (1/simoptions.xTperiod) .* simoptions.drivetimes - (pi/2));
simoptions.Lmode = 0;
simoptions.BuoySim.tether_length = 0;
simoptions.NoOfMachines = 1;

simoptions.ODESim.EvalFcn = 'simplelinearmachineode_proscribedmotion'; 
simoptions.ODESim.PostSimFcn = 'resfun_linear_pscbmot';
simoptions.ODESim.PostPreProcFcn = ['finfun_', mname];

[T, Y, results, design] = simulatemachine_linear(design, simoptions, simoptions.ODESim.PreProcFcn, ...
                                                  simoptions.ODESim.PostPreProcFcn, simoptions.ODESim.EvalFcn, simoptions.ODESim.PostSimFcn); 
                                              
plotresultsproscribedmot_linear(T, Y, results, 1);

%% Test with buoy in sinusoidal sea

% Set up the buoy and sea data files, these are for the 2m buoy
% snappertrunkdir = fileparts(which('wholesystemsim_Snapper'));
% The peak frequency of a PM Specturm
SeaParameters.sigma = 1/3;
% The number of frequencies to be produced
SeaParameters.amp = 0.5;
% The depth of the water
SeaParameters.water_depth = 50;

simoptions.BuoySim.SeaParameters = defaultseaparamaters(SeaParameters);

% use 10 radiation coefficients
simoptions.BuoySim.NRadiationCoefs = 10;

simoptions.BuoySim.tether_length = 6;
simoptions.buoy = 'cyl_4dia_2dr';

simoptions.ODESim.TimeSpan = [0, 60];
% params.amp = 1;
params.sigma = 2 * pi * 0.35;
params.phase = pi/2;

simoptions.BuoySim.SeaParameters = defaultseaparamaters(params);

simoptions.BuoySim.tether_length = 4;

simoptions.ODESim.EvalFcn = 'systemode_linear'; 
simoptions.ODESim.PostPreProcFcn = ['systemfinfun_', mname];
simoptions.ODESim.PostSimFcn = 'systemresfun_linear';    
simoptions.maxAllowedxT = inf;

[T, Y, results, design] = simulatemachine_linear(design, simoptions, simoptions.ODESim.PreProcFcn, ...
                                                  simoptions.ODESim.PostPreProcFcn, simoptions.ODESim.EvalFcn, simoptions.ODESim.PostSimFcn); 

plotresultsbuoysys_linear(T, Y, results, design, 1)


%% Test with buoy in random sea

% Set up the buoy and sea data files, these are for the 2m buoy
% snappertrunkdir = fileparts(which('wholesystemsim_Snapper'));
simoptions.HeaveFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410', 'heave_coefficients_cyl_2di_1dr_d020610.mat');
simoptions.SurgeFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410', 'surge_coefficients_cyl_2di_1dr.mat');
simoptions.HydroCoeffsFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410','cyl_d3103v4.1');
simoptions.ExcitationFile = fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410','cyl_d3103v4.2');
simoptions.BuoySim.BuoyParameters = load(fullfile(getbuoylibdir, 'Cylinder_2m_dia_d010410', 'buoyparams_d3103v4.mat'));
simoptions.buoynum = -1;
design.buoynum = simoptions.buoynum;

simoptions.ODESim.TimeSpan = [0, 60];
% params.amp = 1;
params.peak_freq = 0.35; % centred at resonant frequency
params.phase = pi/2;

simoptions.BuoySim.SeaParameters = defaultseaparamaters(params);

simoptions.BuoySim.tether_length = 4;

simoptions.ODESim.EvalFcn = 'systemode_linear'; 
simoptions.ODESim.PostPreProcFcn = ['systemfinfun_', mname];
simoptions.ODESim.PostSimFcn = 'systemresfun_linear'; 
simoptions.events = 'systemevents_linear';

[T, Y, results, design] = simulatemachine_linear(design, simoptions, simoptions.ODESim.PreProcFcn, ...
                                                  simoptions.ODESim.PostPreProcFcn, simoptions.ODESim.EvalFcn, simoptions.ODESim.PostSimFcn); 

plotresultsbuoysys_linear(T, Y, results, design, 1)