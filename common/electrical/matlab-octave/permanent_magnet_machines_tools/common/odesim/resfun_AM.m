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
% resfun_AM performs a number of post simulation tasks common to all
% machines including the regeneration and post-processing of results.
%
% Input
%
%  T - 
%
%  Y - 
%
%  design - 
%
%  simoptions - structure containing 
%    The following fields must be present in the simoptions structure:
%
%    The following optional fields may also be present in the structure:
%
% See also: oderesults.m
%

    % set some results related options to defaults
    
    simoptions = setfieldifabsent (simoptions, 'CopyPhaseResistance', false);
    % SkipOutputFields: cell array of strings containing fields to be
    % skipped when getting the output values
    simoptions = setfieldifabsent (simoptions, 'SkipOutputFields', {});
    
    % extract the internally calculated results from the ode simulation
    % function. It must be coded so that when called with no arguments it
    % returns a cell array of strings which will become the names of fields
    % in a results structure. The order of the strings must be the order of
    % the return arguments of the ode function when called with more than
    % one input
    results = oderesults(T, Y, simoptions.odeevfun, {design, simoptions}, simoptions.skip, simoptions.SkipOutputFields);
    
    % store the actual simulation time taken
    design.SimTimeSpan = max(T) - simoptions.tspan(1);
    
    % we should use the phase that produced the highest current
    currentInds = simoptions.ODESim.SolutionComponents.PhaseCurrents.SolutionIndices;
    [C,I] = max(max(abs(Y(:,currentInds)), [], 1));
    
    if isfield(results, 'RPhase') &&  simoptions.CopyPhaseResistance
        design.PhaseResistance = results.RPhase;
        results = rmfield(results, 'RPhase');
    else
        results.RPhase = design.PhaseResistance;
    end
    
    % Determine some interesting machine electrical outputs
    design = odeelectricalresults(T, ...
                                  Y(:,currentInds), ...
                                  results.EMF, ...
                                  results.RPhase, ...
                                  design, ...
                                  simoptions);

end
