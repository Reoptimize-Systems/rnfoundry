function delR = Table32r1fDelR(vars)
% function: Table32r1fDelR
% 
% Calculates radial displacement of the circumference of a cylinder
% undergoing a linearly varying radial bod force, as calculated in 'Roark's
% Formulas Stress & Strain 6th edition' in table 32, page 638 case 1 manner
% 1f. 
%
% Input: 
%   
%   vars - (n x 5) matrix:-
%          Col 1. ldeltab, the unit pressure on the vessel (force per unit area)
%          Col 2. a, the outer radius
%          Col 3. b, the inner radius
%          Col 4. E, youngs modulus for the material
%          Col 5. v, Poisson's ratio for the material
%
% Output:
%
%   delR - (n x 2) column vector of values of delR, the radial displacement
%          of the circumference of the cylinder. The first column is the
%          displacement of the outer radius, the second the displacement of
%          the inner radius.
%

    if size(vars,2) == 5
        
        ldeltab = vars(:,1);
        a = vars(:,2);
        b = vars(:,3);
        E = vars(:,4);
        v = vars(:,5);
        
        % first we calculate the displacement of the outer radius (delta a)
        delR(:,1) = ((ldetab .* a) ./ (12 .* E)) .* (((1 - v) .* a.^4) - 8.*(2+v).*a.*b.^3 + 3.*(3 + v).*b.^4 + 6.*(1+v).*a.^2.*b.^2)...
            ./((a-b).*(a.^2 - b.^2));
        
        % delb = max sigma_2 * (b/E), max sigma_2 occurs when r = b
        sigma = ThickWallPressureStresses([ldeltab a b v b], '32.1f');
        
        delR(:,2) =  sigma(:,2) .* (b ./ E);
        
    else
       error('Matrix dimensions do not agree. Table32r1fDelR requires a (n x 5) column matrix') 
    end

end