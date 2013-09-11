function costperkwhr = amortizedenergycostperkwhr(capitalcost, meankw, drate, years, capacityfactor, insurance, oandm)
% calculates the amortized cost of energy over the lifetime of a project
%
% Syntax
%
% costperkwhr = amortizedenergycostperkwhr(capitalcost, meankw, drate, years, capacityfactor)
% costperkwhr = amortizedenergycostperkwhr(..., insurance)
% costperkwhr = amortizedenergycostperkwhr(..., oandm)
%
% Input
%
% capitalcost - capital cost of the project
%
% meankw - mean/rated kW hours produced
% 
% drate - investment discount rate
% 
% years - lifetime of the project in years
% 
% capacityfactor - capacity factor of the device/plant
%
% insurance - annual insurance costs (optional, defaults to zero if not
%   supplied)
%
% oandm - annual operation and maintenace costs (optional, defaults to zero
%   if not supplied)
%
% See also: variableannualcost
%

% Copyright Richard Crozer, The University of Edinburgh

    if nargin < 6
        insurance = 0;
    end
    
    if nargin < 7
        oandm = 0;
    end
    
    % calculare the variable annual cost
    vac = variableannualcost(capitalcost, drate, years);
   
    % calculate the total KW hours produced per year
    totalkwhr = meankw .* capacityfactor .* secondsperyear() ./ (1e3 * 60 * 60);
    
    % finally determine the annual cost of energy (in the capital cost
    % units per kW hour)
    costperkwhr = (vac + insurance + oandm) ./ totalkwhr;
    
end