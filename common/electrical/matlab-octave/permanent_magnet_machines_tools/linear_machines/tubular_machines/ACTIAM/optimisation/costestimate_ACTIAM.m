function [cost, design] = costestimate_ACTIAM(design, simoptions, buoymass)

    if nargin < 3
        buoymass = 0;
    end
    
    [cost, design] = costestimate_TM(design, simoptions, buoymass);
       
end