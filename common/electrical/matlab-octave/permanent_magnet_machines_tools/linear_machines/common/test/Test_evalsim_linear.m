% Test_evalsim_linear

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
design.RlVRp = 10;
design.mu_fT = 0;
design.massT = 0;

design = ratios2dimensions_ACTM(design);

design.PoleWeight = (pi * (design.Rm - design.Rsi)^2 * design.Wp * 7500) + (pi * (design.Ra - design.Ri)^2 * design.Wp * 8600 * design.CoilFillFactor);


simoptions.simfun = 'systemsimfun_ACTM';
mname = 'ACTM';

% Set Parameters

simoptions.Lmode = 0;
simoptions.NoOfMachines = 1;
simoptions.maxAllowedxT = 0.5;

simoptions.buoy = 37;

simoptions.tspan = [0, 60];
% params.amp = 1;

simoptions.tether_length = 4;

simoptions.maxAllowedxT = inf;

simoptions.odeevfun = 'systemode_linear'; 
simoptions.finfun = ['systemfinfun_', mname];
simoptions.resfun = 'systemresfun_linear'; 
simoptions.events = 'systemevents_linear';

params.peak_freq = 1/9; % centred at resonant frequency
params.sigma_range = [0.345575191894877,2.31745966692415;];
params.water_depth = 50;

simoptions.SeaParameters = seasetup('PMPeakFreq', 1/9, ...
                                    'WaterDepth', 50, ...
                                    'NoOfFrequencies', 50);
         
                                
%%

simoptions.DoPreLinSim = true;

[T, Y, results, design, simoptions] = evalsim_linear(design, simoptions);                                


