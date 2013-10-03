function Def = Table3r1dDef(Yvars, E, I, x)
% Table3r1dDef: Calculates the deflection of a beam with both ends fixed,
% as calculated in 'Roark's Formulas for Stress & Strain 6th edition' in
% table 3, row 1d.
%
% Input: 
%   
%   Yvars - (n x 1) column vector of values of R, the radius of the
%          circular cross-section:
%          Yvars(:,1) - W, load at 'a'
%          Yvars(:,2) - l, length of the beam
%          Yvars(:,3) - a, distance from M_A at which 'W' is applied 
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
    
    W = Yvars(:,1);
    l = Yvars(:,2);
    a = Yvars(:,3);
    
    % Calculate the reaction force at A in two stages for readibility
    RA = W .* (l - a).^2 .* (l + 2.*a) ./ (l.^3);
    
    % Calculate the bending moment at A
    MA = -W .* a .* (l - a).^2 ./ (l.^2);
    
    thetaA = 0;
    
    yA = 0;
    
    Def = zeros(size(MA,1),length(x));
    
    for j = 1:size(MA,1)
        for i = 1:length(x)  
            % Calculate the resulting deflection in each case using the
            % generic formula with thetaA = 0 (for fixed supports) and yA =
            % zero (no initial deflection)
            Def(j,i) = GenericYDefConcntdLoad(thetaA, MA(j,1), RA(j,1), W(j,1), x(i), a(j), yA, E, I);
        end
    end
    
end