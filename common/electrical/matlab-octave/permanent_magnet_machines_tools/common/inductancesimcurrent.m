function I = inductancesimcurrent(coilarea, nturns, ampspersqm)
% calculates the current in a coil winding necessary to give a desired
% mean current density over the coil cross-section
%
% Syntax
%
% I = inductancesimcurrent(coilarea, nturns, ampspersqm)
%
% Input
%
% coilarea - cross-sectional area of the coil 
%
% nturns - the number of turns crossing the perpendicular to the
%   cross-section
%
% ampspersqm - the desired mean current density over the coil cross-section
%   in Amps / m^2
%
% Output
%
% I - the current in the coil necessary to give the desired current density

    if nargin < 3
        ampspersqm = 0.1e6;
    end
    
    I = coilarea .* ampspersqm ./ nturns;

end