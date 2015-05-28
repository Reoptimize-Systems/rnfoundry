function [Def, thetaA, MA, RA, yA] = Deflection (Yvars, E, I, x)
% calculates the deflection of a beam with its left end free and its right
% end fixed, undergoing an applied intermediate moment, as calculated in
% 'Roark's Formulas Stress & Strain'.
%
% Input: 
%   
%   Yvars - (n x 1) column vector of values of R, the radius of the
%     circular cross-section:
%          Yvars(:,1) - M0, applied moment at point 'a'
%          Yvars(:,2) - l, length of the beam
%          Yvars(:,3) - a, distance from M_A at which 'M0' is applied 
%
%   E - Young's modulus of the beam material
%
%   I - second moment of inertia of the beam cross-section
%
%   x - row vector of position values at which the deflection is to be
%     calculated
%
% Output:
%
%   Def - (n x 1) column vector of values of the deflection at the
%     corresponding x position
%
    if size (Yvars,2) > 4
        error ('Yvars has too many columns, Yvars must be a (n x 4) matrix')
    end
    
    M0 = Yvars(:,1);
    l = Yvars(:,2);
    a = Yvars(:,3);
    
    % Calculate the angular displacement
    thetaA = -M0 .* (l - a) ./ (E .* I);
    
    % calculate the initial deflection
    yA = M0 .* (realpow (l,2) - realpow (a,2)) ./ (2 .* E .* I);
    
    Def = zeros (size (thetaA,1),length(x));
    
    for j = 1:size(thetaA,1)
        for i = 1:length(x) 
            % Calculate the resulting deflection in each case using the
            % common formula
            % thetaA, MA, RA, M0, x, a, yA, E, I
            Def(j,i) = roark.Beams.ConcMoment.Deflection (thetaA(j,1), 0, 0, M0(j,1), x(i), a(j), yA(j,1), E, I);
        end
    end
    
    if nargout > 1
        RA = zeros (size (thetaA));
        MA = zeros (size (thetaA));
    end

end