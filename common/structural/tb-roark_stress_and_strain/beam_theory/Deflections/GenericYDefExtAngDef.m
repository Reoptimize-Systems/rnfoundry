function yDef = GenericYDefExtAngDef(thetaA, MA, RA, x, a, yA, E, I, theta_0)
% function: GenericYDefExtAngDef
% 
% Calculates the deflection of a beam with an intermediate angular
% deformation using the generic deflection formula for problems of this
% type described in Table 3, row header 4 on page 108 in 'Roark's
% Formulas Stress & Strain 6th edition'. You are required to supply
% information such as the reaction force and moments etc. See Roark for a
% full description of the inputs below.
%
% Input: 
%   
%   thetaA - externally created concentrated angular displacement at point
%            A in radians
%
%   MA - applied couple (moment) at point A
%
%   RA - Reaction force at point A
%
%   x - matrix of position values at which the deflection is to be
%       calculated 
%
%   a - distance from M_A at which angular deformation occurs 
%
%   yA - initial deflection on beam
%
%   E - Young's modulus of the beam material
%
%   I - second moment of inertia of the beam cross-section
%
%   theta_0 - externally created concentrated angular displacement (at 'a')
%
% Output:
%
%   yDef - values of the deflection 
%
    if x < a
        stepfun = 0;
    else
        stepfun = (x-a);
    end
    
    yDef = yA + (thetaA .* x)...
        + ((MA*(x^2))./(2*E*I))...
        + (RA .* (x.^3) ./ (6.*E.*I))...
        + (stepfun.*theta_0);
end