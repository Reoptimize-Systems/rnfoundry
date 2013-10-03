function yDef = GenericYDefConcntdLoad(thetaA, MA, RA, W, x, a, yA, E, I)
% GenericYDefConcntdLoad: Calculates the deflection of a beam with some
% concentrated loading using the generic deflection formula for distributed
% loads described in Table 3, row header row 1 on page 100 in 'Roark's
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
%   W - Load at 'a'
%
%   x - matrix of position values at which the deflection is to be calculated 
%
%   l - length of the beam
%
%   a - distance from M_A at which 'wa' is applied 
%
%   yA - initial deflection on beam
%
%   E - Young's modulus of the beam material
%
%   I - second moment of inertia of the beam cross-section
%
% Output:
%
%   yDef - values of the deflection 
%
    if x < a
        stepfun = 0;
    else
        stepfun = (x-a).^3;
    end
    
    yDef = yA + (thetaA .* x)...
        + ((MA*(x^2))./(2*E*I))...
        + (RA .* (x.^3) ./ (6.*E.*I))...
        - (stepfun.*W./(6.*E.*I));
end