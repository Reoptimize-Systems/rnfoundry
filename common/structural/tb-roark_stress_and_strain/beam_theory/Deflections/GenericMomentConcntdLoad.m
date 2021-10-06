function M = GenericMomentConcntdLoad(MA, RA, W, x, a)
% GenericMomentConcntdLoad: Calculates the moments in a beam with some
% concentrated loading using the generic deflection formula for distributed
% loads described in Table 3, row header row 1 on page 100 in 'Roark's
% Formulas Stress & Strain 6th edition'. You are required to supply
% information such as the reaction force and moments etc. See Roark for a
% full description of the inputs below.
%
% Input: 
%   
%   MA - applied couple (moment) at point A
%
%   RA - Reaction force at point A
%
%   W - Load at 'a'
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
    stepfun = zeros (size(x));
    
    stepfun(x>=a) = x(x>=a) - a(x>=a);
    
    M = MA + RA .* x - W .* stepfun;
    
end