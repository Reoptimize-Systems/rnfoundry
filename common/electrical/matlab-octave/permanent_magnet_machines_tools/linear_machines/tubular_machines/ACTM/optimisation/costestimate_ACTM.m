function [cost, design] = costestimate_ACTM(design, simoptions, buoymass)
% extimates the cost of the air-cored tubular linear machine
%
% Syntax
%
% [cost, design] = costestimate_ACTM(design, simoptions, buoymass)
% 

    if nargin < 3
        buoymass = 0;
    end
    
    [cost, design] = costestimate_TM(design, simoptions, buoymass);
    
end