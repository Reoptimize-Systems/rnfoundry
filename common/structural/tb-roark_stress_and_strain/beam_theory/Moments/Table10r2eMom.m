function Mom = Table10r2eMom(Yvars, E, I, x)
% Table10r2eMom: Calculates the moments in a beam with its left end simply
% supported and its right end simply supported, undergoing a linearly
% distributed load, and an axial load at one end as calculated in 'Roark's
% Formulas Stress & Strain 6th edition' in table 10, page 166 row 2e.
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
%   Mom - (n x 1) column vector of values of the deflection at the
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
    
    k = sqrt(P./(E.*I));
    
    % Calculate the reaction force at A
    RA = Table10r2eRA(wa, wl, l, a);
    
    % Calculate the angular displacement in two stages
    thetaA = Table10r2eThetaA(P, wa, wl, l, a, k);
    
    Mom = zeros(size(thetaA,1),length(x));
    
    for j = 1:size(thetaA,1)
            % Calculate the resulting moments in each case using the
            % generic formula with MA = 0 (for simply supported) 
            Mom(j,:) = GenericMomentDistribTransLoadAndAxialLoad(thetaA(j,1), 0, RA(j,1), P(j,1), wa(j,1), wl(j,1), E, I, l(j,1), a(j,1), x);
    end
    
end