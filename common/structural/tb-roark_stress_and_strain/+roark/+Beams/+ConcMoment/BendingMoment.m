function M = BendingMoment (MA, RA, M0, x, a)
% Calculates the moments in a beam with an applied concentrated
% intermediate moment using the generic deflection formula for concentrated
% loads described in 'Roark's Formulas Stress & Strain'. You are required
% to supply information such as the reaction force and moments etc. See
% Roark for a full description of the inputs below.
%
% Input: 
%   
%   MA - applied couple (moment) at point A
%
%   RA - Reaction force at point A
%
%   M0 - Applied moment at 'a'
%
%   x - matrix of position values at which the deflection is to be calculated 
%
%   l - length of the beam
%
%   a - distance from M_A at which 'W' is applied 
%
% Output:
%
%   M - moment 
%

    stepfun = zeros (size (x));
    
    stepfun(x>=a) = 1;
    
    M = MA + RA .* x + M0 .* stepfun;
    
end