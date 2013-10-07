e%% Test_designandevaluate_PMSM

% First we need some plausible machine variables for testing, we will use
% the AWS variables for this

Phases = 3;         % Number of Phases in machine
Taup = 0.1;         % Pole Pitch
bp = 0.8*Taup;      % Magnet Pitch
lm = 15/1000;       % Magnet depth
kfill = 0.585;      % Copper fill factor
g = 3/1000;         % Initial air-gap
dc = 1/1000;        % 1 mm diameter wire for both machines
E = [200e9 151e9];  % Young's modulus of elasticity for Structural steel and laminated steel respectively
hty = 30/1000;      % Translator/Field yoke (back iron) thickness
hsy = hty;          % Stator/Armature yoke (half central section) thickness
Taus = Taup/Phases; % Combined tooth and slot height
bs = Taus/2;        % Slot height
bt = Taus/2;        % Tooth height
hs = 100/1000;      % Tooth depth
options.targetPower = 5e5; % 500kW machine
options.mlength = 6; % Overlap between stator and translator, i.e. stator is mleng metres longer than the translator
options.pointsPerPole = 30;
%%
% First do a stack length of 1.5 m
ls = 1.5; 
% Calculate the extra length needed for the bearings
bearingWidth = 0.1; 
options.alphab = (ls + 2*bearingWidth) / ls;

% Get the dimensionless ratios from the parameters
[dgVlm, dtVdg, bpVlm, taupVbp, lsVbp, btVbc, dcVbs, dtiVdt, dbiVlm] = dimensions2ratios_PMSM(lm, Taup, bp, g, hs, ls, bt, dc, hty, lm+g+hs+hsy);

RgVRc = 10; % Ratio of machine resistance to grid resistance

[out, design1] = designandevaluate_PMSM(dgVlm, dtVdg, bpVlm, taupVbp, lsVbp, btVbc, dcVbs, dtiVdt, dbiVlm, RgVRc, kfill, lm, options);

%%
% Then do a stack length of 0.5 m
ls = 0.5; 
bearingWidth = 0.1; 
options.alphab = (ls + 2*bearingWidth) / ls;

% Get the dimensionless ratios from the parameters
[dgVlm, dtVdg, bpVlm, taupVbp, lsVbp, btVbc, dcVbs, dtiVdt, dbiVlm] = dimensions2ratios_PMSM(lm, Taup, bp, g, hs, ls, bt, dc, hty, lm+g+hs+hsy);

RgVRc = 10; % Ratio of machine resistance to grid resistance

[out, design2] = designandevaluate_PMSM(dgVlm, dtVdg, bpVlm, taupVbp, lsVbp, btVbc, dcVbs, dtiVdt, dbiVlm, RgVRc, kfill, lm, options);

%% Original Design

% First do a stack length of 1 m
ls = 1; 
% Calculate the extra length needed for the bearings
bearingWidth = 0.1; 
options.alphab = (ls + 2*bearingWidth) / ls;

% Get the dimensionless ratios from the parameters
[dgVlm, dtVdg, bpVlm, taupVbp, lsVbp, btVbc, dcVbs, dtiVdt, dbiVlm] = dimensions2ratios_PMSM(lm, Taup, bp, g, hs, ls, bt, dc, hty, lm+g+hs+hsy);

RgVRc = 10; % Ratio of machine resistance to grid resistance

[out, design3] = designandevaluate_PMSM(dgVlm, dtVdg, bpVlm, taupVbp, lsVbp, btVbc, dcVbs, dtiVdt, dbiVlm, RgVRc, kfill, lm, options);

%% Minimum Stack length
ls = 0.1; 
bearingWidth = 0.1; 
options.alphab = (ls + 2*bearingWidth) / ls;

% Get the dimensionless ratios from the parameters
[dgVlm, dtVdg, bpVlm, taupVbp, lsVbp, btVbc, dcVbs, dtiVdt, dbiVlm] = dimensions2ratios_PMSM(lm, Taup, bp, g, hs, ls, bt, dc, hty, lm+g+hs+hsy);

RgVRc = 10; % Ratio of machine resistance to grid resistance

[out, design_min_stack] = designandevaluate_PMSM(dgVlm, dtVdg, bpVlm, taupVbp, lsVbp, btVbc, dcVbs, dtiVdt, dbiVlm, RgVRc, kfill, lm, options);