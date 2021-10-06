function Def = Deflection (Yvars, E, I, x)
% Calculates the deflection of a beam with both ends simply supported, as
% calculated in 'Roark's Formulas Stress & Strain'.
%
% Input: 
%   
%   Yvars - (n x 1) column vector of values of R, the radius of the
%     circular cross-section:
%          Yvars(:,1) - W, load at 'a'
%          Yvars(:,2) - l, length of the beam
%          Yvars(:,3) - a, distance from M_A at which 'W' is applied 
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
    if size(Yvars,2) > 3
        error('Yvars has too many columns, Yvars must be a (n x 3) matrix')
    end
    
    W = Yvars(:,1);
    l = Yvars(:,2);
    a = Yvars(:,3);
    
    thetaA = roark.Beams.ConcLoad.LURF.thetaA (W,l,a,E,I);
    
    yA = roark.Beams.ConcLoad.LURF.yA (W,l,a,E,I);
    
    Def = zeros(size(W,1),length(x));
    
    for j = 1:size(W,1)
        for i = 1:length(x)  
            % Calculate the resulting deflection in each case using the
            % generic formula with thetaA = 0 (for fixed supports) and yA =
            % zero (no initial deflection)
            Def(j,i) = roark.Beams.ConcLoad.Deflection (thetaA(j,1), 0, 0, W(j,1), x(i), a(j), yA(j,1), E, I);
        end
    end
    
end