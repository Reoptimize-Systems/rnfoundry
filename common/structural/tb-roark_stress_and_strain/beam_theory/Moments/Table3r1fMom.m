function Mom = Table3r1fMom(Yvars, x)
% Table3r1fDef: Calculates the deflection of a beam with its left end
% guided and its right end simply supported, as calculated in 'Roark's
% Formulas Stress & Strain 6th edition' in table 3, page 101 row 1f.
%
% Input: 
%   
%   Yvars - (n x 1) column vector of values of R, the radius of the
%          circular cross-section:
%          Yvars(:,1) - W, load at 'a'
%          Yvars(:,2) - l, length of the beam
%          Yvars(:,3) - a, distance from M_A at which 'W' is applied 
%
%   x - row vector of position values at which the deflection is to be calculated 
%
% Output:
%
%   Mom - (n x 1) column vector of values of the moment at the
%         corresponding x position
%
    if size(Yvars,2) > 3
        error('Yvars has too many columns, Yvars must be a (n x 4) matrix')
    end
    
    W = Yvars(:,1);
    l = Yvars(:,2);
    a = Yvars(:,3);
    
    % Reaction force at A is zero
    RA = 0;
    
    % Calculate the bending moment at A
    MA = Table3r1fMA(W, l, a);
    
    Mom = zeros(size(MA,1),length(x));
    
    for j = 1:size(MA,1)
        for i = 1:length(x)  
            % Calculate the resulting moment in each case using the
            % generic formula
            Mom(j,i) = GenericMomentConcntdLoad(MA(j), RA(j), W(j), x(i), a(j));
        end
    end
    
end