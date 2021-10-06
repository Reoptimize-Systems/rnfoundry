function Def = Table10r2eDef(Yvars, E, I, x)
% Calculates the deflection of a beam with its left end simply supported
% and its right end simply supported, undergoing a linearly distributed
% load, and an axial load at one end as calculated in 'Roark's Formulas
% Stress & Strain 6th edition' in table 10, page 166 row 2e.
%
% Input: 
%   
%   Yvars - (n x 1) column vector of values of R, the radius of the
%          circular cross-section:
%          Yvars(:,1) - P, axial load at end A
%          Yvars(:,2) - wa, unit load at 'a'
%          Yvars(:,3) - wl, unit load at M_B, the end of the beam
%          Yvars(:,4) - l, length of the beam
%          Yvars(:,5) - a, distance from M_A at which 'wa' is applied 
%
%   E - Young's modulus of the beam material
%
%   I - second moment of inertia of the beam cross-section
%
%   x - row vector of position values at which the deflection is to be calculated 
%
% Output:
%
%   Def - (n x 1) column vector of values of the deflection at the
%   corresponding x position
%
    if size(Yvars,2) > 5
        error('Yvars has too many columns, Yvars must be a (n x 5) matrix')
    end
    
    P = abs(Yvars(:,1));
    wa = Yvars(:,2);
    wl = Yvars(:,3);
    l = Yvars(:,4);
    a = Yvars(:,5);
    
    k = sqrt(P ./ (E.*I));
    
    % Calculate the reaction force at A
    RA = Table10r2eRA(wa, wl, l, a);
    
    % Calculate the angular displacement
    thetaA = Table10r2eThetaA(P, wa, wl, l, a, k);
    
    Def = zeros(size(thetaA,1),length(x));
    
    yA = 0; 
    MA = 0;
    
    for j = 1:size(thetaA,1)
            % Calculate the resulting deflection in each case using the
            % generic formula with MA = 0 (for simply supported) and yA =
            % zero (no initial deflection)
            Def(j,:) = GenericYDefDistribTransLoadAndAxialLoad(thetaA(j,1), yA, MA, RA(j,1), P(j,1), wa(j,1), wl(j,1), E, I, l(j,1), a(j,1), x);
            
    end
    
end