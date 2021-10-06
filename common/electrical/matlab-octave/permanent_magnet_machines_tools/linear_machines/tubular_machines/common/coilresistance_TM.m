function design = coilresistance_TM (design)
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

    if ~isfield (design, 'MTL')
        
        % for a tubular machine the MTL is the radial distance to the coil
        % 2D shape centroid
        if isfield (design, 'CoilCentroid')
            design.MTL = design.CoilCentroid(1);
        else
            % if not avaialble use the 
            design.MTL = 2 * pi * design.Rcm;
        end

    end
    
    design.CoilResistance = wireresistancedc ('round', design.Dc, design.MTL*design.CoilTurns);

end