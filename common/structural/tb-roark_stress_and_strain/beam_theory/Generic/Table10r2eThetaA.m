function thetaA = Table10r2eThetaA(P, wa, wl, l, a, k)
% Table10r2eRA: Calculates the slope at end A in a beam with its left end
% simply supported and its right end simply supported, undergoing a
% linearly distributed load, and an axial load at one end as calculated in
% 'Roark's Formulas Stress & Strain 6th edition' in table 10, page 166 row
% 2e.
%
% Input: 
%
%   P - axial load at end A
%
%   wa - unit load at 'a'
%
%   wl - unit load at M_B, the end of the beam
%
%   l - length of the beam
%
%   a - distance from M_A at which 'wa' is applied 
%
%   k - the value of the function (P/EI)^(0.5)
%
% Output:
%
%   thetaA - (n x 1) column vector of values of the slope thetaA at end A
%
    thetaA = (-wa ./ (k .* P)) .* ...
                ( ( (1 - cos(k .* (l - a))) ./ sin(k .* l) ) - ...
                (( k .* ( (l - a).^2 )) ./ (2 .* l)) );
    
    thetaA = thetaA - ( (wl - wa) ./ (k .* P) ) .* ...
                ( ((k .* (l - a) - sin( k .* (l - a) )) ./ (k .* (l - a) .* sin(k .* l))) - ...
                (k .* ((l - a).^2) ./ (6 .* l)) );
    
end