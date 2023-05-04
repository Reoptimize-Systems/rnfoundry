function [results, design, summary_time_inds] = resfun_RADIAL_SLOTTED(T, Y, design, simoptions)
% post processes results from an ode simulation of a slotted radial flux
% electrical machine
%
% Syntax
%
% [results, design] = resfun_RADIAL_SLOTTED(T, Y, design, simoptions)
%
%
% Input
%
% 

% Copyright Richard Crozier 2014


    [results, design, summary_time_inds] = resfun_RADIAL(T, Y, design, simoptions);
    
    if isfield(results, 'TqaddEBD') && (size (results.TqaddEBD, 2) > 4)
        
        % get iron losses forces
        results.TqCogging = results.TqaddEBD(:,5);
        
    end
    
end