function MA = MA (W, l, a)
% Calculates the Moment at end A of a beam with its left end guided and its
% right end simply supported, as calculated in 'Roark's Formulas Stress &
% Strain'
%
% Input: 
%   
%   W - load at 'a'
%
%   l - length of the beam
%
%   a - distance from M_A at which 'W' is applied 
%
% Output:
%
%   MA - (n x 1) column vector of values of the deflection at the
%     corresponding x position
%

    % Calculate the bending moment at A
    MA = W.*(l - a);
    
end