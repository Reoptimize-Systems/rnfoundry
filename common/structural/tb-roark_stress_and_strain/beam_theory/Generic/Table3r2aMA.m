function MA = Table3r2aMA()
% Table3r2dMA: Calculates the moment at end A for a beam with its left end
% fixed and its right end fixed, as calculated in 'Roark's Formulas Stress
% & Strain 6th edition' in table 3, page 104 row 2d.
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
%   Mom - values of MA
%
    % Calculate the bending moment at A in two stages
    MA = 0;
    
end