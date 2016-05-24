function [T, Y, results, design, simoptions] = simulatemachine_linear (design, simoptions, varargin)
% performs a simulation of a linear machine design and system operation
% using ode solvers
%
% Syntax
%
% [T, Y, results, design, simoptions] = simulatemachine_linear (design, ...
%                   simoptions, 'Parameter', 'Value')
%
% Inputs
%
% design, simoptions - structures containing all the information necessary
% to perform the machine simultion. Their contents will depend on the
% particular simulation and design being used
%

% Copyright Richard Crozer, The University of Edinburgh

    % call the more generic function simulatemachine_AM
    [T, Y, results, design, simoptions] = ...
        simulatemachine_AM (design, simoptions, varargin{:});

end