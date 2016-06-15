%% ACTM

% This section sets up an air cored tubular machine design in preparation
% for simulation
clear design simoptions

% Number of phases in machine
design.Phases = 3;  
% Magnet outer radius
design.Rm = 0.1;
% air gap
design.g = 5/1000;
% Coil inner radius
design.Ri = design.Rm + design.g;
% ratio of magnet width to pole width
design.WmVWp = 0.75;
% ratio of pole width to magnet outer radius
design.WpVRm = 0.5;
% ratio of coil inner radius to magnet outer radius
design.RiVRm = design.Ri / design.Rm;
% ratio of coil outer radius to magnet outer radius
design.RoVRm = 1.2;
% ratio of coil sheath outer radius to coil outer radius
design.RaVRo = 1.025;
% ratio of magnet stack shaft outer radius to magnet outer radius
design.RsoVRm = 0.1;
% ratio of magnet stack shaft inner radius to magnet stack shaft outer
% radius
design.RsiVRso = 0;
% ratio of coil axial height to pole height
design.WcVWp = 1/3;
% the coil wire configuration
design.CoilFillFactor = 0.65;
design.CoilTurns = 50;

design.mode = 2; 

% poles on the translator and stator respectively
design.Poles = [10 30];
% series/parallel branch configuration in phase winding
design.CoilsPerBranch = 10;
design.Branches = 1;

% convert the ratio based description to actual physical dimensions
design = ratios2dimensions_ACTM(design);


%% Set up Common Parameters

% ratio of load resistance to phase resistance
design.RlVRp = 10;
% ratio of load inductance to phase inductance
design.LgVLc = 0;

%% Test with buoy in 1m amplitude sinusoidal sea

% Set up the buoy and sea data files, these are for the 2m buoy with 1m
% draft. There are two ways to specify a buoy, in this case we will use a
% string, which refers to a buoy directory name in the library of buoys
simoptions.buoy = 'cyl_2dia_1dr';

% create a sea, in this case a single frequency sea of frequency 0.35 Hz
% and default amplitude 0.5m
simoptions.SeaParameters = seasetup ('Sigmas', 2 * pi * 0.35);

% the system is simulated as a heaving buoy attached to a tether which
% passes through a hole in a vertically mounted structure and raises the
% translator vertically. The tether_length is the initial distance from the
% tether point on the buoy (the centre of the buoy) to the hole when the
% buoy is at rest in the starting position.
simoptions.tether_length = 4;

% set the end stop position (the maximum allowed translator displacement).
% Here we set it to inf, so there are no end stops and the translator can
% travel any distance
simoptions.maxAllowedxT = inf;

% choose a time span for the simulation
simoptions.tspan = [0, 60];

% set up some functions to be used to generate the data and run the
% simulation
simoptions.simfun = 'systemsimfun_ACTM';
simoptions.odeevfun = 'systemode_linear'; 
simoptions.finfun = 'systemfinfun_ACTM';
simoptions.resfun = 'systemresfun_linear'; 
simoptions.events = 'systemevents_linear'; 

% run the simulation with the specified parameters
[T, Y, results, outdesign, outsimoptions] = simulatemachine_AM(design, simoptions); 

% create several plots showing different outputs from the simulation
plotresultsbuoysys_linear(T, Y, results, design, outsimoptions, 1)


%% SImulate with buoy in random sea

% Set up the buoy and sea data files, here we specify a buoy using the
% alternative method, an integer. This selects the 37th buoy in the library
% directory
simoptions.buoy = 37;

simoptions.tether_length = 4;

simoptions.maxAllowedxT = inf;

simoptions.simfun = 'systemsimfun_ACTM';
simoptions.odeevfun = 'systemode_linear'; 
simoptions.finfun = 'systemfinfun_ACTM';
simoptions.resfun = 'systemresfun_linear'; 
simoptions.events = 'systemevents_linear';

% this time we specify a pierson-moskowitz sea state with a peak frequency
% of 1/9
simoptions.SeaParameters = seasetup ('PMPeakFreq', 1/9);

% we'll do a 60 second simulation, but for a 
simoptions.tspan = [0, 60];

[T, Y, results, outdesign, outsimoptions] = simulatemachine_linear (design, simoptions); 
    
% produce some interesting plots of the output, but plotting only every 5th
% value in the output arrays for efficiency
skip = 5;
plotresultsbuoysys_linear(T, Y, results, outdesign, outsimoptions, skip)


 
 
 