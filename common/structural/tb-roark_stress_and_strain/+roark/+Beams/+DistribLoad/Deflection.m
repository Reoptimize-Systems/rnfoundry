function yDef = Deflection (thetaA, MA, RA, wa, wl, x, l, a, yA, E, I)
% Calculates the deflection of a beam with some distributed loading using
% the generic deflection formula for distributed loads described in
% 'Roark's Formulas Stress & Strain'. You are required to supply
% information such as the reaction force and moments etc. See Roark for a
% full description of the inputs below.
%
% Input: 
%   
%   thetaA - externally created concentrated angular displacement at end
%     A in radians
%
%   MA - applied couple (moment) at point A
%
%   RA - Reaction force at point A
%
%   wa - unit load at 'a'
%
%   wl - unit load at M_B, the end of the beam
%
%   x - matrix of position values at which the deflection is to be
%     calculated
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
        stepfun1 = 0;
        stepfun2 = 0;
    else
        stepfun1 = realpow (x-a,4);
        stepfun2 = realpow (x-a,5);
    end
    
    yDef = yA + (thetaA .* x)...
        + ((MA*(x^2))/(2*E*I))...
        + (RA .* (x.^3) ./ (6.*E.*I))...
        - (stepfun1.*wa./(24.*E.*I))...
        -(stepfun2 .* (wl-wa) ./ (120.*E.*I.*(l-a)));
end