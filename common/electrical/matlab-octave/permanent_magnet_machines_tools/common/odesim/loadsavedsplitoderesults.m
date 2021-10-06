function [T, Y, results] = loadsavedsplitoderesults(basename, rootdir)
% loads results saved to disk at each split point during a split ODE
% simulation 
%
% Syntax
%
% [T, Y, results] = loadsavedsplitoderesults(basename, rootdir)
%
% Descriptions
%
% loadsavedsplitoderesults loads a set of structures saved to .mat files
% which contain results gathered at the intermediate points of a simulation
% using the matlab ode solvers with the odesplit.m function. The names of
% the files containing the results must have the format:
%
% basename_#.mat
%
% where basename is a starting string chosen by the user and # represents
% an integer which counts upwards from 1 at each split point. Each file
% should contain a solution structure of the same format as produced by the
% ode solver routines, and a second structure containing other internally
% produced results of interest named odeinternals. 
%
% The x and y fields of the solution structure will be concatenated to
% produce T and Y vectors for the full duration of the simulation. Each
% field of the odeinternals structure will be be concatenated to produce
% the results structure for the full simulation.
%
%

    % rootdir = homedir;
    if nargin < 2
        rootdir = pwd;
    end

    T = [];
    Y = [];

    thefiles = dir(fullfile(rootdir, [basename, '_*.mat']));

    if isempty(thefiles)
        error('No split output files found in directory.')
    end
    
    for find = 1:numel(thefiles)

        load(fullfile(rootdir, [basename, '_', num2str(find), '.mat']), 'odeinternals', 'sol');

        if find == 1
            
            T = [T; sol.x'];
            
            Y = [Y; sol.y'];
        
            results = odeinternals;
            fnames = fieldnames(results);
        else
            
            T = [T; sol.x(2:end)'];
            
            Y = [Y; sol.y(:,2:end)'];
            
            for ii = 1:numel(fnames)
                results.(fnames{ii}) = [results.(fnames{ii}); odeinternals.(fnames{ii})(2:end,:)];
            end
        end

    end  
    
end