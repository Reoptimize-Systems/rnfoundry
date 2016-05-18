%% Test_split_ode_sim

close all
clear design NSdesign Sdesign NSsimoptions Ssimoptions S_NS_res_design

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

%%

NSsimoptions.NoOfMachines = 2;

% set up the functions for a no split ode evaluation
NSsimoptions.simfun = 'simfun_ACTIAM';
NSsimoptions.finfun = 'systemfinfun_ACTIAM'; 
NSsimoptions.odeevfun = 'systemode_linear'; 
NSsimoptions.resfun = 'systemresfun_linear';

% use buoy number 37, 4m diameter, 2m draft
NSsimoptions.buoy = 'cyl_4dia_2dr';

% NSsimoptions.SeaParameters = seasetup('sigma', 2 * pi * 0.35, ...
%                                     'phase', pi / 2);

NSsimoptions.SeaParameters = seasetup('PMPeakFreq', 1/9);

NSsimoptions.tether_length = 5;
NSsimoptions.water_depth = 40;

% simoptions.ODESim.InitialConditions = [0,0,0];
NSsimoptions.skip = 1;
NSsimoptions.tspan = [0, 90];
% simoptions.Lmode = 0;

%% Run split 

% run the simulation
[NST, NSY, NSresults, NSdesign, NSsimoptions] = simulatemachine_linear(design, NSsimoptions, ...
                                                             NSsimoptions.simfun, ...
                                                             NSsimoptions.finfun, ... 
                                                             NSsimoptions.odeevfun, ...
                                                             NSsimoptions.resfun);

%% run split system

Ssimoptions.NoOfMachines = NSsimoptions.NoOfMachines;

% set up the functions for a no split ode evaluation
Ssimoptions.simfun = NSsimoptions.simfun;
Ssimoptions.finfun = NSsimoptions.finfun;
Ssimoptions.odeevfun = NSsimoptions.odeevfun;
Ssimoptions.resfun = 'splitsystemresfun_linear';  
Ssimoptions.spfcn = 'splitodesystemres_linear';
Ssimoptions.splitode = 10;
Ssimoptions.SaveSplitResults = true;

% use buoy number 37, 4m diameter, 2m draft
Ssimoptions.buoy = NSsimoptions.buoy;

Ssimoptions.SeaParameters = NSsimoptions.SeaParameters;
Ssimoptions.tether_length = NSsimoptions.tether_length;
Ssimoptions.water_depth = NSsimoptions.water_depth;

Ssimoptions.tspan = NSsimoptions.tspan;

%% Run split 

% run the simulation
[ST, SY, Sresults, Sdesign, Ssimoptions] = simulatemachine_linear(design, Ssimoptions, ...
                                                             Ssimoptions.simfun, ...
                                                             Ssimoptions.finfun, ... 
                                                             Ssimoptions.odeevfun, ...
                                                             Ssimoptions.resfun);
                                                         
% load the saved split time series results from disk
[ST, SY, STresults] = loadsavedsplitoderesults(char(Ssimoptions.spfcn));

% post process the results from the ode sim using the non split version of
% the results function, and the results loaded from disk
[S_NS_res_results, S_NS_res_design] = systemresfun_linear(ST, SY, Sdesign, Ssimoptions);

%% Compare

% rowheadings = {'PowerLoadMean', ...
%                'PowerLoadPeak', ...
%                'EnergyLoadTotal', ...
%                'PowerPhaseRMean', ...
%                'PowerSystemMean', ...
%                'EMFPhaseRms', ...
%                'IPhaseRms', ...
%                'JCoilRms', ...
%                'IPhasePeak', ...
%                'JCoilPeak', ...
%                'vRmax'};

rowheadings = unique([fieldnames(Sdesign); fieldnames(NSdesign); fieldnames(S_NS_res_design)]);

data = [];
rmrowheadings = [];
for i = 1:numel(rowheadings)
    if isnumeric(NSdesign.(rowheadings{i}))
        data = [data; NSdesign.(rowheadings{i})(1), ...
                      Sdesign.(rowheadings{i})(1), ...
                      S_NS_res_design.(rowheadings{i})(1), ...
                      100 * abs((NSdesign.(rowheadings{i})(1) - Sdesign.(rowheadings{i})(1))/NSdesign.(rowheadings{i})(1))];
    else
        rmrowheadings = [rmrowheadings, i];
    end
end

% remove nonumeric rowheadings
rowheadings(rmrowheadings) = [];

colheadings = {'Non Split', 'Split', 'Split Recalc', 'S-NS diff %'};

displaytable(data, colheadings, 15, [], rowheadings);

plotresultsbuoysys_linear(NST, NSY, NSresults, NSdesign);
% 
% plotresultsbuoysys_linear(ST, SY, STresults, Sdesign);

