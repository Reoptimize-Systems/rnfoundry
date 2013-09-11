function [dc, fulldc] = ConductorDiameter(coilArea, ksfill, Ntot)
% ConductorDiameter: a function to calculate the diameter of wire necessary
% to achieve a desired number of turns in a coil of a given cross-sectional
% area and fill factor. 
%
% Syntax
%
% [dc, fulldc] = ConductorDiameter(coilArea, ksfill, Ntot)
%
% Description
%
% The conductor area necessary takes account of the coating or sheath on
% the wire
% Input:
%   
%   coilArea - the total cross-sectional area of the coil in m^2
%
%   kfill - the fill factor of the coil
%
%   Ntot = the desired number of turns in the coil
%
% Output:
%
%   dc - the conductor diameter (in m) necessary to achieve the turn number
%   specified in Ntot
%
%   fulldc - the total diameter of the wire including sheath
%
    coilArea = coilArea * ksfill;
    
    % calculate necessary conductor diameter in mm to achieve turns
    maxdc = 1000 * sqrt((4 * coilArea ./ (Ntot * pi)));
    
    dc = maxdc;
    
    % we will reduce dc until the conductor + sheath/enamel is less than or
    % equal to the necessary maximum to achieve the turns
    fulldc = maxdc * 1.01;
    
    while fulldc > maxdc
        
        % reduce dc by 1%
        dc = dc * 0.99;

        % determine the full cross-sectional diameter
        if dc < 1.6
            % First polynomial covers conductor diameter from 0.1 mm to 1.6 mm.
            fulldc = dc * (1 + (-0.1135*dc^5 + 1.6055*dc^4 - 8.5416*dc^3 + 21.481*dc^2 - 27.039*dc + 18.561)/100);
        else
            % For greater than 1.6 mm, a power fit was used
            fulldc = dc * (1 + (5.9131 * dc^(-0.6559))/100); 
        end
       
    end
    
    % convert back to m
    dc = dc / 1000;
    
    fulldc = fulldc / 1000;
    
end