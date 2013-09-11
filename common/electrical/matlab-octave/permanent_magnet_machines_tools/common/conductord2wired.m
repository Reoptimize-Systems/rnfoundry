function fulldc = conductord2wired(conductordc)
% converts conductor diameter to a total wire diamter including the sheath
% or enamel
%
% Syntax
%
% fulldc = conductord2wired(conductordc)
%
% Description
%
% Converts the conductor diameter, specified in metres to the total wire
% diameter that this would result in in metres when the outer sheath, or
% enamel coating is included
%
% see data from FD sims for more information on how the diameter is
% calculated
%
%

    % Convert dc to mm for polynomials (poly input variables are in mm)
    dc = conductordc * 1000;
    
    if dc < 1.6
        % First polynomial covers conductor diameter from 0.1 mm to 1.6 mm.
        fulldc = dc * (1 + (-0.1135*dc^5 + 1.6055*dc^4 - 8.5416*dc^3 + 21.481*dc^2 - 27.039*dc + 18.561)/100);
    else
        % For greater than 1.6 mm, a power fit was used so it can be
        % extrapolated to larger diameter wire 
        fulldc = dc * (1 + (5.9131 * dc^(-0.6559))/100); 
    end
    
    % convert the conductor diameter from mm to m
    fulldc = fulldc / 1000;
    
end