function thetaA = Table3r2aThetaA(wa, wl, l ,a, E, I)
% Table3r1fRA: Calculates the initial deflection at end A of a beam with
% its left end guided and its right end simply supported, as calculated in
% 'Roark's Formulas Stress & Strain 6th edition' in table 3, page 101 row
% 1f.
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
%   E - Young's modulus of the beam material
%
%   I - second moment of inertia of the beam cross-section
%
% Output:
%
%   thetaA - (n x 1) column vector of values of the deflection
%
    
    thetaA = (wa .* (l - a).^3 ./ (6.*E.*I)) + ((wl - wa) .* (l - a).^3 ./ (24.*E.*I));
    
end