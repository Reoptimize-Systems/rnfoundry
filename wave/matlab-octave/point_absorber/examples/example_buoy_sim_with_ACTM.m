%% ACTM

% This section sets up an air cored tubular machine design in preparation
% for simulation
clear design simoptions

design = test_design_TM_SLOTLESS ();

%% Set up Common Parameters

% ratio of load resistance to phase resistance
design.RlVRp = 10;
% ratio of load inductance to phase inductance
design.LgVLc = 0;

simoptions.MagFEASimType = 'multiple';

%% Test with buoy in 1m amplitude sinusoidal sea

% Set up the buoy and sea data files, these are for the 2m buoy with 1m
% draft. There are two ways to specify a buoy, in this case we will use a
% string, which refers to a buoy directory name in the library of buoys
simoptions.BuoySim.buoy = 'cyl_2dia_1dr';

% create a sea, in this case a single frequency sea of frequency 0.35 Hz
% and default amplitude 0.5m
simoptions.BuoySim.SeaParameters = seasetup ('Sigmas', 2 * pi * 0.35, ...
                                             'Amplitudes', 1);

% the system is simulated as a heaving buoy attached to a tether which
% passes through a hole in a vertically mounted structure and raises the
% translator vertically. The tether_length is the initial distance from the
% tether point on the buoy (the centre of the buoy) to the hole when the
% buoy is at rest in the starting position.
simoptions.BuoySim.tether_length = 4;

% set the end stop position (the maximum allowed translator displacement).
% Here we set it to inf, so there are no end stops and the translator can
% travel any distance
simoptions.BuoySim.maxAllowedxT = inf;

% choose a time span for the simulation
simoptions.ODESim.TimeSpan = [0, 60];

% set up some functions to be used to generate the data and run the
% simulation
simoptions.ODESim.PreProcFcn = 'simfun_TM_SLOTLESS';
simoptions.ODESim.EvalFcn = 'systemode_linear'; 
simoptions.ODESim.PostPreProcFcn = 'systemfinfun_TM_SLOTLESS';
simoptions.ODESim.PostSimFcn = 'systemresfun_linear'; 
simoptions.events = 'systemevents_linear'; 

design = completedesign_TM_SLOTLESS (design, simoptions);

% run the simulation with the specified parameters
[T, Y, results, outdesign, outsimoptions] = simulatemachine_AM(design, simoptions); 

% create several plots showing different outputs from the simulation
plotresultsbuoysys_linear(T, Y, results, design, outsimoptions, 1)


%% SImulate with buoy in random sea

% Set up the buoy and sea data files, here we specify a buoy using the
% alternative method, an integer. This selects the 37th buoy in the library
% directory
simoptions.buoy = 37;

simoptions.BuoySim.tether_length = 4;

simoptions.maxAllowedxT = inf;

simoptions.ODESim.PreProcFcn = 'systemsimfun_ACTM';
simoptions.ODESim.EvalFcn = 'systemode_linear'; 
simoptions.ODESim.PostPreProcFcn = 'systemfinfun_ACTM';
simoptions.ODESim.PostSimFcn = 'systemresfun_linear'; 
simoptions.events = 'systemevents_linear';

% this time we specify a pierson-moskowitz sea state with a peak frequency
% of 1/9
simoptions.BuoySim.SeaParameters = seasetup ('PMPeakFreq', 1/9);

% we'll do a 60 second simulation, but for a 
simoptions.ODESim.TimeSpan = [0, 60];

[T, Y, results, outdesign, outsimoptions] = simulatemachine_linear (design, simoptions); 
    
% produce some interesting plots of the output, but plotting only every 5th
% value in the output arrays for efficiency
skip = 5;
plotresultsbuoysys_linear(T, Y, results, outdesign, outsimoptions, skip)


 
 
 