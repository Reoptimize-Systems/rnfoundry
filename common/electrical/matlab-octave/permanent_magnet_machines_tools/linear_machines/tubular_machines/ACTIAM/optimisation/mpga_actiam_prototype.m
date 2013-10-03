% mpga_actm_prototype.m         (Multi Population Genetic Algorithm)
%
% This script implements the Multi Population Genetic Algorithm.
% Real valued representation for the individuals is used.
%
% The purpose of the optimisation is to find the lowest cost per watt of an
% air-cored tubular machine for the Industrial Electroics Journal paper.
%
% Author:     Richard Crozier
% History:    05.02.2009     file created
% -------------------------------------------------------------------------
% Input Variables:

GGAP = .8;           % Generation gap, how many new individuals are created
INSR = .9;           % Insertion rate, how many of the offspring are inserted
XOVR =  1;           % Crossover rate
SP = 2;              % Selective Pressure
MUTR = 1;            % Mutation rate; only a factor;
MIGR = 0.2;          % Migration rate between subpopulations
MIGGEN = 20;         % Number of generations between migration (isolation time)

TERMEXACT = 1e-4;    % Value for termination if minimum reached

SEL_F = 'sus';       % Name of selection function
XOV_F = 'recint';    % Name of recombination function for individuals
MUT_F = 'mutbga';    % Name of mutation function
OBJ_F = 'objactiam_prototype';   % Name of function for objective values

% -------------------------------------------------------------------------
%
% Get boundaries of objective function
FieldDR = feval(OBJ_F,[],1);

% compute SUBPOP, NIND depending on number of variables (defined in objective function)
NVAR = size(FieldDR,2);           % Get number of variables from objective function
SUBPOP = 2;                       % Number of subpopulations
NIND = 50;                        % Number of individuals per subpopulations
MAXGEN = 180; % Max number of generations
MUTR = MUTR / NVAR;           % Mutation rate depending on NVAR

STEP = 1;
DISPLAYMODE = 2;
SAVEMODE = 1;

%
% [Best, IndAll, Chrom, ObjV] = mpgafun(OBJ_F, GGAP, INSR, XOVR, SP, MIGR, MIGGEN, TERMEXACT, SEL_F,...
%     XOV_F, MUT_F, SUBPOP, NIND, MAXGEN, MUTR, STEP, DISPLAYMODE, SAVEMODE)

filename = fullfile(fundir('mpga_actiam_prototype'), 'objactiam_prototype_output_4.mat');

resumeFile = 'objactiam_prototype_output_4.mat';

[Best, IndAll, Chrom, ObjV] = mpgafun(OBJ_F, GGAP, INSR, XOVR, SP, MIGR, MIGGEN, TERMEXACT, SEL_F,...
    XOV_F, MUT_F, SUBPOP, NIND, MAXGEN, MUTR, STEP, DISPLAYMODE, SAVEMODE, filename, resumeFile);




