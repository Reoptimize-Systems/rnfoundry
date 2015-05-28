function M = BendingMoment (MA, RA, M0, x, a)
% Calculates the moments in a beam with an applied concentrated
% intermediate moment using the generic described in 'Roark's Formulas
% Stress & Strain'. You are required to supply information such as the
% reaction force and moments etc. See Roark for a full description of the
% inputs below.
%
% Input: 
%   
%   MA - applied couple (moment) at end A
%
%   RA - Reaction force at end A
%
%   M0 - Applied moment at 'a'
%
%   x - matrix of position values at which the deflection is to be
%     calculated
%
%   a - distance from M_A at which 'M0' is applied 
%
% Output:
%
%   M - moment at each point in x
%

    stepfun = zeros (size (x));
    
    stepfun(x>=a) = 1;
    
    M = MA + RA .* x + M0 .* stepfun;
    
end