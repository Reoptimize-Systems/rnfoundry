function Slope = Table3r1fSlope(Yvars, E, I, x)
% function: Table3r1fSlope
% 
% Calculates the deflection of a beam with its left end guided
% and its right end simply supported, as calculated in 'Roark's Formulas
% Stress & Strain 6th edition' in table 3, page 101 row 1f.
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
    RA = 0;
    
    % Calculate the bending moment at A
    MA = W.*(l - a);
    
    thetaA = 0;
    
    Slope = zeros(size(MA,1),length(x));
    
    for j = 1:size(MA,1)
        for i = 1:length(x)  
            % Calculate the resulting slope in each case using the
            % generic formula with thetaA = 0 (for fixed supports) and RA =
            % zero 
            Slope(j,i) = GenericSlopeConcntdLoad(thetaA, MA(j,1), RA, W(j,1), x(i), a(j,1), E, I);
        end
    end
    
end

