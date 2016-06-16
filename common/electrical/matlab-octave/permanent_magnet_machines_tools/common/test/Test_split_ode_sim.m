%% Test_split_ode_sim

close all
clear design NSdesign Sdesign NSsimoptions Ssimoptions S_NS_res_design

design = test_design_ACTM ();

design = ratios2dimensions_ACTIAM(design);

%%

NSsimoptions.NoOfMachines = 2;

% set up the functions for a no split ode evaluation
NSsimoptions.ODESim.PreProcFcn = 'simfun_ACTIAM';
NSsimoptions.ODESim.PostPreProcFcn = 'systemfinfun_ACTIAM'; 
NSsimoptions.ODESim.EvalFcn = 'systemode_linear'; 
NSsimoptions.ODESim.PostSimFcn = 'systemresfun_linear';

% use buoy number 37, 4m diameter, 2m draft
NSsimoptions.BuoySim.buoy = 'cyl_4dia_2dr';

% NSsimoptions.BuoySim.SeaParameters = seasetup('sigma', 2 * pi * 0.35, ...
%                                     'phase', pi / 2);

NSsimoptions.BuoySim.SeaParameters = seasetup('PMPeakFreq', 1/9);

NSsimoptions.BuoySim.tether_length = 5;

% simoptions.ODESim.InitialConditions = [0,0,0];
NSsimoptions.skip = 1;
NSsimoptions.ODESim.TimeSpan = [0, 90];
% simoptions.Lmode = 0;

%% Run split 

% run the simulation
[NST, NSY, NSresults, NSdesign, NSsimoptions] = simulatemachine_linear(design, NSsimoptions);

%% run split system

Ssimoptions.NoOfMachines = NSsimoptions.NoOfMachines;

% set up the functions for a no split ode evaluation
Ssimoptions.ODESim.PreProcFcn = NSsimoptions.ODESim.PreProcFcn;
Ssimoptions.ODESim.PostPreProcFcn = NSsimoptions.ODESim.PostPreProcFcn;
Ssimoptions.ODESim.EvalFcn = NSsimoptions.ODESim.EvalFcn;
Ssimoptions.ODESim.PostSimFcn = 'splitsystemresfun_linear';  
Ssimoptions.ODESim.SplitPointFcn = 'splitodesystemres_linear';
Ssimoptions.ODESim.Split = 10;
Ssimoptions.ODESim.SaveSplitResults = true;

% use buoy number 37, 4m diameter, 2m draft
Ssimoptions.BuoySim.buoy = NSsimoptions.BuoySim.buoy;

Ssimoptions.BuoySim.SeaParameters = NSsimoptions.BuoySim.SeaParameters;
Ssimoptions.BuoySim.tether_length = NSsimoptions.BuoySim.tether_length;

Ssimoptions.ODESim.TimeSpan = NSsimoptions.ODESim.TimeSpan;

%% Run split 

% run the simulation
[ST, SY, Sresults, Sdesign, Ssimoptions] = simulatemachine_linear(design, Ssimoptions);
                                                         
% load the saved split time series results from disk
[ST, SY, STresults] = loadsavedsplitoderesults(char(Ssimoptions.ODESim.SplitPointFcn));

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
    if isfield (NSdesign, rowheadings{i}) && isfield (Sdesign, rowheadings{i}) && isfield (S_NS_res_design, rowheadings{i})
        if isnumeric(NSdesign.(rowheadings{i}))
            data = [data; NSdesign.(rowheadings{i})(1), ...
                          Sdesign.(rowheadings{i})(1), ...
                          S_NS_res_design.(rowheadings{i})(1), ...
                          100 * abs((NSdesign.(rowheadings{i})(1) - Sdesign.(rowheadings{i})(1))/NSdesign.(rowheadings{i})(1))];
        else
            rmrowheadings = [rmrowheadings, i];
        end
    else
        rmrowheadings = [rmrowheadings, i];
    end
end

% remove nonumeric rowheadings
rowheadings(rmrowheadings) = [];

colheadings = {'Non Split', 'Split', 'Split Recalc', 'S-NS diff %'};

displaytable(data, colheadings, 15, [], rowheadings);

% T, Y, results, design, simoptions, skip
plotresultsbuoysys_linear(NST, NSY, NSresults, NSdesign, NSsimoptions);
% 
% plotresultsbuoysys_linear(ST, SY, STresults, Sdesign);

