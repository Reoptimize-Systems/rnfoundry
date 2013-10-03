function Def = Table3r2dDef(Yvars, E, I, x)
% Table3r2dDef: Calculates the deflection of a beam with its left end fixed
% and its right end fixed, as calculated in 'Roark's Formulas Stress &
% Strain 6th edition' in table 3, page 104 row 2d.
%
% Input: 
%   
%   Yvars - (n x 1) column vector of values of R, the radius of the
%          circular cross-section:
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
%   corresponding x position
%
    if size(Yvars,2) > 4
        error('Yvars has too many columns, Yvars must be a (n x 4) matrix')
    end
    
    wa = Yvars(:,1);
    wl = Yvars(:,2);
    l = Yvars(:,3);
    a = Yvars(:,4);
    
    % Calculate the reaction force at A in two stages for readibility
    RA = (wa ./ ((l.^3).*2)) .* ((l-a).^3) .* (l+a);
    
    RA = RA + (((wl-wa) ./ (20 .* (l.^3)))  .* ((l-a).^3) .* ((3.*l) + (2.*a)) );
    
    % Calculate the bending moment at A in two stages
    MA = (-wa ./ ((l.^2).*12)) .* ((l-a).^3) .* (l+(3.*a));
    
    MA = MA - (((wl-wa) ./ (60 .* (l.^2))) .* ((l-a).^3) .* ((2.*l)+(3.*a)));
    
    Def = zeros(size(MA,1),length(x));
    
    for j = 1:size(MA,1)
        for i = 1:length(x)  
            % Calculate the resulting deflection in each case using the
            % generic formula with thetaA = 0 (for fixed supports) and yA =
            % zero (no initial deflection)
            Def(j,i) = GenericYDefDistribLoad(0, MA(j,1), RA(j,1), wa(j,1), wl(j,1), x(i), l(j), a(j), 0, E, I);
        end
    end
    
end