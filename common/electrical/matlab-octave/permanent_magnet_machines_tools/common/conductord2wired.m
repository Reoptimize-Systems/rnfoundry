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
    
    fulldc = ones (size (dc)) * nan;
    
    % First polynomial covers conductor diameter from 0.1 mm to 1.6 mm.
    fulldc(dc < 1.6) = dc(dc < 1.6) ...
        .* ( 1 + (-0.1135.*dc(dc < 1.6).^5 ...
                 + 1.6055.*dc(dc < 1.6).^4 ...
                 - 8.5416.*dc(dc < 1.6).^3 ...
                 + 21.481.*dc(dc < 1.6).^2 ...
                 - 27.039.*dc(dc < 1.6) ...
                 + 18.561 ) ...
               ./ 100);
    
    % For greater than 1.6 mm, a power fit was used so it can be
    % extrapolated to larger diameter wire 
    fulldc(dc >= 1.6) = dc(dc >= 1.6) .* (1 + (5.9131 * dc(dc >= 1.6).^(-0.6559))./100); 

    % convert the conductor diameter from mm to m
    fulldc = fulldc ./ 1000;
    
end