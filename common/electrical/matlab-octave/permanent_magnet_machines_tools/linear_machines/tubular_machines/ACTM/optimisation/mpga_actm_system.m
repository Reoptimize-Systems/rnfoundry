% mpga_actm_prototype.m         (Multi Population Genetic Algorithm)
%
% This script implements the Multi Population Genetic Algorithm.
% Real valued representation for the individuals is used.
%
% The purpose of the optimisation is to find the lowest cost per watt of an
% air-cored tubular machine for the Industrial Electroics Journal paper.
%
% Author:     Richard Crozier
% -------------------------------------------------------------------------
% Input Variables:

common_setup_options_opt_and_comp_linear_gens_for_wave_energy

% pathstr = fileparts(which('common_setup_options_opt_and_comp_linear_gens_for_wave_energy'));

% load(fullfile(pathstr, 'setup_options_opt_and_comp_linear_gens.mat'));

mpgaoptions.OBJ_F = 'objactm_system';   % Name of function for objective values

% set up the functions
simoptions.simfun = 'dummysimfun_ACTM';
simoptions.finfun = 'systemfinfun_ACTM';
simoptions.odeevfun = 'systemode_linear';
% simoptions.dpsidxfun = 'polypsidot_ACTIAM'; %@dpsidx_tubular; 
simoptions.resfun = 'splitsystemresfun_linear';
simoptions.splitode = 10;
simoptions.spfcn = 'splitodesystemres_linear';
simoptions.events = 'systemevents_linear';

simoptions.DisplayDesignFcn = 'displaydesign_ACTM';

%% Set up the design evaluation options

options.E = [207e9 100e5 207e9];
%options.targetPower = 10e3; % 10kW machine
%options.mlength = 4; % Overlap between stator and translator, i.e. stator is mleng metres longer than the translator
% options.pointsPerPole = 40;
options.coilYieldStrength = 70e6;

%% Perform the optimisation

ObjectiveArgs = {simoptions, options};

% Get boundaries of objective function
FieldDR = feval(mpgaoptions.OBJ_F,[],1, simoptions, options);

% compute SUBPOP, NIND depending on number of variables (defined in objective function)
mpgaoptions.NVAR = size(FieldDR,2);           % Get number of variables from objective function
mpgaoptions.SUBPOP = 3;                       % Number of subpopulations
mpgaoptions.NIND = 30;                        % Number of individuals per subpopulations
mpgaoptions.MAXGEN = 100; % Max number of generations
mpgaoptions.MUTR = mpgaoptions.MUTR / mpgaoptions.NVAR;  % Mutation rate depending on NVAR

mpgaoptions.STEP = 1;
mpgaoptions.DISPLAYMODE = 2;
mpgaoptions.SAVEMODE = 1;

%
% [Best, IndAll, Chrom, ObjV] = mpgafun(OBJ_F, GGAP, INSR, XOVR, SP, MIGR, MIGGEN, TERMEXACT, SEL_F,...
%     XOV_F, MUT_F, SUBPOP, NIND, MAXGEN, MUTR, STEP, DISPLAYMODE, SAVEMODE)

mpgaoptions.filename = 'objactm_system_output.mat';
% 
mpgaoptions.RESUME = true;
mpgaoptions.resumefile = mpgaoptions.filename;

% [Best, IndAll, Chrom, ObjV] = mpgafun(OBJ_F, GGAP, INSR, XOVR, SP, MIGR, MIGGEN, TERMEXACT, SEL_F,...
%     XOV_F, MUT_F, SUBPOP, NIND, MAXGEN, MUTR, STEP, DISPLAYMODE, SAVEMODE, filename)

[Best, IndAll, Chrom, ObjV, mpgaoptions] = mpgafun2(mpgaoptions, 'ObjectiveArgs', ObjectiveArgs);

%% Evaluate the winner

% set up the functions
simoptions.simfun = 'simfun_ACTM';
simoptions.finfun = 'systemfinfun_ACTM';
simoptions.odeevfun = 'systemode_linear';
simoptions.resfun = 'systemresfun_linear';

if isfield(simoptions, 'splitode')
    simoptions = rmfield(simoptions, 'splitode');
end

if isfield(simoptions, 'spfcn')
    simoptions = rmfield(simoptions, 'spfcn');
end

simoptions.events = 'systemevents_linear';

% [design, simoptions] = preprocsystemdesign_ACTM(simoptions, BestInd(1:14));

[design, simoptions] = preprocsystemdesign_ACTM(simoptions, IndAll(end,:));

[score, design, simoptions, T, Y, results] = designandevaluate_ACTM(design, simoptions, options);

plotresultsbuoysys_linear(T, Y, results, design, 1)

