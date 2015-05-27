function [Def, thetaA, MA, RA, yA] = Deflection (Yvars, E, I, x)
% Calculates the deflection of a beam with its left end simply supported
% and its right end simply supported, undergoing a linearly distributed
% load, as calculated in 'Roark's Formulas Stress & Strain'.
%
% Input: 
%   
%   Yvars - (n x 1) column vector of values of R, the radius of the
%     circular cross-section:
%
%          Yvars(:,1) - wa, unit load at 'a'
%          Yvars(:,2) - wl, unit load at M_B, the end of the beam
%          Yvars(:,3) - l, length of the beam
%          Yvars(:,4) - a, distance from M_A at which 'wa' is applied 
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
%     corresponding x position
%
    if size (Yvars,2) > 4
        error ('Yvars has too many columns, Yvars must be a (n x 4) matrix')
    end
    
    wa = Yvars(:,1);
    wl = Yvars(:,2);
    l = Yvars(:,3);
    a = Yvars(:,4);
    
    % Calculate the reaction force at A in two stages for readibility
    RA = (wa ./ (l.*2)) .* ((l-a).^2);
    
    RA = RA + (((l-a).^2) .* (wl-wa) ./ (6 .* l));
    
    % Calculate the angular displacement in two stages
    thetaA = (-wa ./ (24 .* E .* l.* I)) .* ((l-a).^2) .* (l.^2 + (2.*l.*a) - a.^2);
    
    thetaA = thetaA - ((wl-wa) ./ (360 .* E .* l.* I)).*((l-a).^2) .*( (7.*(l.^2)) + (6.*l.*a) - (3.*(a.^2)));
    
    Def = zeros (size (thetaA,1),length(x));
    
    for j = 1:size(thetaA,1)
        for i = 1:length(x) 
            % Calculate the resulting deflection in each case using the
            % generic formula with MA = 0 (for simply supported) and yA =
            % zero (no initial deflection)
            Def(j,i) = roark.Beams.DistribLoad.Deflection (thetaA(j,1), 0, RA(j,1), wa(j,1), wl(j,1), x(i), l(j), a(j), 0, E, I);
        end
    end
    
    if nargout > 1
        MA = zeros (size (thetaA));
        yA = zeros (size (thetaA));
    end

end