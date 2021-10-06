% MPGA_IEPO.m         (Multi Population Genetic Algorithm)
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

options.GGAP = .8;           % Generation gap, how many new individuals are created
options.INSR = .9;           % Insertion rate, how many of the offspring are inserted
options.XOVR =  1;           % Crossover rate

options.SEL_F = 'sus';       % Name of selection function
options.XOV_F = 'recint';    % Name of recombination function for individuals
options.MUT_F = 'mutbga';    % Name of mutation function
options.OBJ_F = 'objactm_fixedspeed';   % Name of function for objective values

% -------------------------------------------------------------------------
%
% Get boundaries of objective function
FieldDR = feval(options.OBJ_F,[],1);

% compute SUBPOP, NIND depending on number of variables (defined in objective function)
options.NVAR = size(FieldDR,2);           % Get number of variables from objective function
options.SUBPOP = 4;                       % Number of subpopulations
options.NIND = 25;                        % Number of individuals per subpopulations
options.MAXGEN = 200;                     % Max number of generations

options.STEP = 5;
options.DISPLAYMODE = 2;
options.SAVEMODE = 1;
options.filename = 'objactm_fixedspeed_output.mat';

%options.resumeFile = options.filename;

[Best, IndAll, Chrom, ObjV, options] = mpgafun2(options);


                                  
                                  
                                  