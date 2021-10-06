function RA = Table10r2eRA(wa, wl, l, a)
% Table10r2eRA: Calculates the reation force at end A in a beam with its
% left end simply supported and its right end simply supported, undergoing
% a linearly distributed load, and an axial load at one end as calculated
% in 'Roark's Formulas Stress & Strain 6th edition' in table 10, page 166
% row 2e.
%
% Input: 
%   
%   wa - unit load at 'a'
%
%   wl - unit load at M_B, the end of the beam
%
%   l - length of the beam
%
%   a - distance from M_A at which 'wa' is applied 
%
% Output:
%
%   RA - (n x 1) column vector of values of the reaction force RA at end A
%
    RA = (wa ./ (l.*2)) .* ((l - a).^2);
    
    RA = RA + (((l-a).^2) .* (wl-wa) ./ (6 .* l));
    
end