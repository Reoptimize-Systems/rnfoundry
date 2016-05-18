%% Test_simulatemachine_ACTIAM

clear design simoptions

design.Phases = 3;         % Number of Phases in machine
design.Rm = 0.15;
design.g = 3/1000;
design.Ri = design.Rm + design.g;
design.WmVWp = 0.75;
design.WpVRm = 0.4;
design.RiVRm = design.Ri / design.Rm;
design.RoVRm = 1.2;
design.RaVRo = 1.03;
design.RsoVRm = 0.1;
design.RsiVRso = 0;
design.WcVWp = 1/3;
design.Rs2VHmag = 0.5;
design.Rs1VHmag = 0.5;
design.Ws2VhalfWs = 0.5;
design.Ws1VhalfWs = 0.5;

design.CoilFillFactor = 0.55;
design.Dc = 0.5/1000;        % 1 mm diameter wire 
design.mode = 2;
design.LgVLc = 0;
design.RlVRp = 10; % Ratio of machine resistance to grid resistance

design.Poles = [5 10];

design = ratios2dimensions_ACTIAM(design);

% set up the functions
simoptions.simfun = 'simfun_ACTIAM';
simoptions.finfun = 'systemfinfun_ACTIAM'; 
simoptions.odeevfun = 'systemode_linear'; 
simoptions.resfun = 'systemresfun_linear';

% use buoy number 37, 4m diameter, 2m draft
simoptions.buoy = 'cyl_4dia_2dr';

simoptions.SeaParameters = seasetup('sigma', 2 * pi * 0.35, ...
                                    'phase', pi / 2);
simoptions.tether_length = 5;
simoptions.water_depth = 40;

% simoptions.ODESim.InitialConditions = [0,0,0];
simoptions.skip = 1;
simoptions.tspan = [0, 10];
% simoptions.Lmode = 0;

%%

% run the simulation
[T, Y, results, design, simoptions] = simulatemachine_linear(design, simoptions, ...
                                                            simoptions.simfun, ...
                                                            simoptions.finfun, ... 
                                                            simoptions.odeevfun, ...
                                                            simoptions.resfun);

%%
% plot the results
plotresultsbuoysys_linear(T, Y, results, design);
