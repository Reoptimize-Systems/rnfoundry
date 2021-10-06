function Def = Table10r2fDef(Yvars, E, I, x)
% function: Table10r2fDef
% 
% Calculates the deflection of a beam with its left end guided
% and its right end simply supported, undergoing a linearly distributed
% load, and an axial load at one end as calculated in 'Roark's Formulas
% Stress & Strain 6th edition' in table 10, page 166 row 2f.
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
    
    P = Yvars(:,1);
    wa = Yvars(:,2);
    wl = Yvars(:,3);
    l = Yvars(:,4);
    a = Yvars(:,5);
    
    [Cn, k] = AxialLoadCCoeffs(P', E', I', l');
    
    Cn = Cn';
    
    k = k';
    
    Can = AxialLoadCaCoeffs(P', E', I', a', l');
    
    Can = Can';
    
    % Calculate the reaction force at A in two stages for readibility
    MA = wa.*Can(:,3)./(k.^2.*Cn(:,1)) + (wl - wa).*Can(:,4)./(k.^3.*(l - a).*Cn(:,1));
    
    % Calculate the angular displacement in two stages
    yA = -wa.*(Can(:,3)./Cn(:,1) - k.^2.*(l-a).^2./2)./(k.^2.*P);
    
    Def = zeros(size(yA,1),length(x));
    
    for j = 1:size(yA,1)
            % Calculate the resulting deflection in each case using the
            % generic formula with MA = 0 (for simply supported) and yA =
            % zero (no initial deflection)
            Def(j,:) = GenericYDefDistribTransLoadAndAxialLoad(0, yA(j,1), MA(j,1), 0, P(j,1), wa(j,1), wl(j,1), E, I, l(j,1), a(j,1), x);
    end
    
end