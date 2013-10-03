function R = coilresistance_TM(RoVRm, Rm, g, N, rwire)
% Calculates the resistance of a coil in a tubular machine
%
% Arguments: (input)
%
%   RoVRm - Ro /Rm ratio
%
%   Rm - The translator radius, can be in any metric units provided the
%   value of the air-gap is given in the same units.
%
%   g - The length of the air-gap, must be in the same units as Rm.
%
%   N - Number of turns in coil
%
%   rwire - radius of wire used
%
%  Arguments: (output)
%
%   R - The total resistance of a single coil
%
% 
% See also: wirelength_TM.m

% Created by Richard Crozier 2013

    length = wirelength_TM(RoVRm, N, Rm, g);
    
    R = 1.68e-8 * length / (pi*rwire^2);

end