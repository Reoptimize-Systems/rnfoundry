function RB = Table3r2aRB(wa, wl, l ,a)
% Table3r2dMA: Calculates the reaction force at end A for a beam with its
% left end fixed and its right end fixed, as calculated in 'Roark's
% Formulas Stress & Strain 6th edition' in table 3, page 104 row 2d.
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
%   RB - values of RB
%

    % Calculate the reaction force at A 
    RB = (wa + wl) .* (l - a) ./ 2;
end