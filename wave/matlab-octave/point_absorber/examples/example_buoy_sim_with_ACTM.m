%% ACTM

clear design simoptions

design.Phases = 3;         % Number of phases in machine
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
design.CoilTurns = 50;
design.mode = 2; 
design.LgVLc = 0;
design.Poles = [10 30];
design.CoilsPerBranch = 10;
design.Branches = 1;
% design.FieldDirection = 1;
% design.PowerPoles = poles(1);

design = ratios2dimensions_ACTM(design);


%% Set up Common Parameters

design.RlVRp = 10;

simoptions.Lmode = 0;
simoptions.NoOfMachines = 1;
simoptions.maxAllowedxT = 0.5;

%% Test with buoy in sinusoidal sea

% Set up the buoy and sea data files, these are for the 2m buoy
simoptions.buoy = 'cyl_2dia_1dr';

simoptions.tspan = [0, 60];
% params.amp = 1;
% params.sigma = 2 * pi * 0.35;
% params.phase = pi/2;

% use a default sea
simoptions.SeaParameters = seasetup ('Sigmas', 2 * pi * 0.35);

% set the end stop position (the maximum allowed translator displacement)
simoptions.maxAllowedxT = inf;

simoptions.tether_length = 4;

simoptions.simfun = 'systemsimfun_ACTM';
simoptions.odeevfun = 'systemode_linear'; 
simoptions.finfun = 'systemfinfun_ACTM';
simoptions.resfun = 'systemresfun_linear'; 
simoptions.events = 'systemevents_linear'; 

[T, Y, results, outdesign, outsimoptions] = simulatemachine_linear(design, simoptions, simoptions.simfun, ...
                                                  simoptions.finfun, simoptions.odeevfun, simoptions.resfun); 

plotresultsbuoysys_linear(T, Y, results, design, outsimoptions, 1)


%% SImulate with buoy in random sea

% Set up the buoy and sea data files, these are for the 2m buoy
simoptions.buoy = 37;
% design.buoynum = simoptions.buoynum;

simoptions.tspan = [0, 60];
% params.amp = 1;

simoptions.tether_length = 4;

simoptions.maxAllowedxT = inf;

simoptions.simfun = 'systemsimfun_ACTM';
simoptions.odeevfun = 'systemode_linear'; 
simoptions.finfun = 'systemfinfun_ACTM';
simoptions.resfun = 'systemresfun_linear'; 
simoptions.events = 'systemevents_linear';

% params.peak_freq = 1/9; % centred at resonant frequency
% params.sigma_range = [0.345575191894877,2.31745966692415;];
% params.water_depth = 50;

simoptions.SeaParameters = seasetup ('PMPeakFreq', 1/9);

[T, Y, results, outdesign, outsimoptions] = simulatemachine_linear (design, simoptions); 
    
% produce some interesting plots of the output, but plotting only every 5th
% value in the output arrays for efficiency
skip = 5;
plotresultsbuoysys_linear(T, Y, results, outdesign, outsimoptions, skip)


 
 
 