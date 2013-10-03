function Def = Table3r4dDef(Yvars, E, I, x)
% function: Table3r4dDef
% 
% Calculates the deflection of a beam with its left end fixed and its right
% end fixed, undergoing an externally created angular deformation, as
% calculated in 'Roark's Formulas Stress & Strain 6th edition' in table 3,
% page 108 row 4d.
%
% Input: 
%   
%   Yvars - (n x 1) column vector of values of R, the radius of the
%          circular cross-section:
%          Yvars(:,1) - theta_0, externally created angular displacement at 'a'
%          Yvars(:,2) - l, length of the beam
%          Yvars(:,3) - a, distance from M_A at which theta_0 is applied 
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
    if size(Yvars,2) > 3
        error('Yvars has too many columns, Yvars must be a (n x 4) matrix')
    end
    
    theta_0 = Yvars(:,1);
    l = Yvars(:,2);
    a = Yvars(:,3);
    
    % Calculate the reaction force at A in two stages for readibility
    RA = 6.*E.*I.*theta_0.*(l-2.*a) ./ (l.^3);
    
    % Calculate the bending moment at A
    MA = 2.*E.*I.*theta_0.*(3.*a - 2.*l) ./ (l.^2);
    
    Def = zeros(size(MA,1),length(x));
    
    for j = 1:size(MA,1)
        for i = 1:length(x)  
            % Calculate the resulting deflection in each case using the
            % generic formula with thetaA = 0 (for fixed supports) and yA =
            % zero (no initial deflection)
            Def(j,i) = GenericYDefExtAngDef(0, MA(j,1), RA(j,1), x(i), a(j,1), 0, E, I, theta_0(j,1));
        end
    end
    
end