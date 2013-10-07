function [results, design] = resfun_AM(T, Y, design, simoptions)
% calculates the results from a generator/system ode simulation performed
% using the appropriatately coded simulation functions
%
% Syntax
%
% [results, design] = resfun_AM(T, Y, design, simoptions)
%
% Description
%
% TODO: add description
%
%
% See also: oderesults.m
%

    % extract the internally calculated results from the ode simulation
    % function. It must be coded so that when called with no arguments it
    % returns a cell array of strings which will become the names of fields
    % in a results structure. The order of the strings must be the order of
    % the return arguments of the ode function when called with more than
    % one input
    results = oderesults(T, Y, simoptions.odeevfun, {design, simoptions}, simoptions.skip);
    
    % store the actual simulation time taken
    design.SimTimeSpan = max(T) - simoptions.tspan(1);
    
    % we should use the phase that produced the highest current
    [C,I] = max(max(abs(Y(:,simoptions.ODEPhaseCurrentCol:(simoptions.ODEPhaseCurrentCol-1+design.Phases))), [], 1));
    
    if isfield(results, 'RPhase')
        design.PhaseResistance = results.RPhase;
        results = rmfield(results, 'RPhase');
    end
    
    % Determine some interesting machine electrical outputs
    design = odeelectricalresults(T, ...
                                  Y(:,simoptions.ODEPhaseCurrentCol:simoptions.ODEPhaseCurrentCol-1+design.Phases), ...
                                  results.EMF, ...
                                  design, ...
                                  simoptions);

end
