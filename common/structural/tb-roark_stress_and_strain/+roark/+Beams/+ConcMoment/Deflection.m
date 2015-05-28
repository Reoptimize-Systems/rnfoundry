function yDef = Deflection (thetaA, MA, RA, M0, x, a, yA, E, I)
% calculates the deflection of a beam with an applied concentrated
% intermediate moment based on the formula presented in 'Roark's Formulas
% Stress & Strain'. You are required to supply information such
% as the reaction force and moments etc. See Roark for a full description
% of the inputs below.
%
% Input: 
%   
%   thetaA - externally created concentrated angular displacement at end
%     A in radians
%
%   MA - applied couple (moment) at end A
%
%   RA - Reaction force at end A
%
%   M0 - applied couple at a
%
%   x - matrix of position values at which the deflection is to be calculated 
%
%   a - distance from M_A at which 'M0' is applied 
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
        stepfun = realpow (x-a,2);
    end
    
    yDef = yA + (thetaA .* x) ...
        + ( MA .* realpow (x,2) ./ (2.*E.*I) ) ...
        + ( RA .* realpow (x,3) ./ (6.*E.*I) ) ...
        + ( M0 .* stepfun ./ (2.*E.*I) );
end