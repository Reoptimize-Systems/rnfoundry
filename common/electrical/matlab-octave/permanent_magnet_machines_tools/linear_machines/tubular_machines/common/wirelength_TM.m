function length = wirelength_TM(RoVRm, N, Rm, g)
% Calculates the total length of wire in a coil in a tubular machine in the
% same units as the input dimensions.
%
% Arguments: (input)
%
%   RoVRm - Ro / Rm ratio
%
%   N - Number of turns in coil
%
%   Rm - The translator radius, can be in any metric units provided the
%        value of the air-gap is given in the same units.
%
%   g - The length of the air-gap, must be in the same units as Rm.
%
%  Arguments: (output)
%
%   length - total length of wire in one coil of machine, this must be
%            multiplied by the total number of coils in phase to get
%            complete length of wire in phase.
%
%
% See also: coilresistance_TM

% Created by Richard Crozier 2013


    %coil height
    ch = (RoVRm * Rm) - g - Rm;
    
    % Average Height assuming even distribution
    ch = ch / 2;
    
    % Average coil radius
    ch = Rm + g + ch;
    
    % Average turn length
    length = pi * 2 * ch;
    
    % Total coil length
    length = length * N;

end