function [N, dc] = CoilTurns(area, CoilFillFactor, dc)
% calculates the number of turns in an available area of coil
% cross-section.
%
% Syntax
%
% [N, dc] = CoilTurns(area, CoilFillFactor, dc)
%
% Input:
%   
%   area - the total cross-sectional area of the coil in m^2
%
%   CoilFillFactor - the fill factor of the coil
%
%   dc - the target conductor diameter in m, this may be modified to the
%        closest that is actually possible
%
% Output:
%
%   N - the possible number of turns achievable for that area
%
%   dc - the conductor diameter achieved for this number of turns
%

    % First determine the full wire thickness including the insulation, for
    % details on these polynomials see page 124 of lab book 2. The
    % polynomials actually return the percentage increase in wire thickness
    % due to the insulation.
    
    fulldc = conductord2wired(dc);

    % Calculate the full conductor cross-sectional area, including the
    % insulation
    Ac = pi * (fulldc/2)^2;
    
    % determine if we can fit any turns in
    if Ac > (area * CoilFillFactor)
        % If we cannot, we must modify the conductor diameter until at
        % least one does. 
        N = 1;
        % Allow 1% for insulation, this calculation of dc will return a
        % value in m, so no need to convert back
        dc = 0.99 * sqrt(4 * area * CoilFillFactor / pi);
    else
        % If we can fit turns in, calculate how many
        N = round((area * CoilFillFactor) / Ac);
    end
    
end