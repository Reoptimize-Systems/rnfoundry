function delR = Table32r1eDelR(vars)
% function: Table32r1eDelR
% 
% Calculates radial displacement of the circumference of a cylinder
% undergoing a uniformally distributed radial body force, as calculated in 'Roark's
% Formulas Stress & Strain 6th edition' in table 32, page 638 case 1 manner
% 1e. 
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
        delR(:,1) = (ldeltab .* a.^3 ./ (3 .* E)) .* (1 - v + 2 .* (2 + v) .* b.^2 ./ (a .* (a + b)));
        
        delR(:,2) =  (ldeltab .* a .* b ./ (3 .* E)) .* ((b ./ a) .* (1 - v) + 2 .* a .* (2 + v) ./ (a + b));
        
    else
       error('Matrix dimensions do not agree. Table32r1eDelR requires a (n x 5) column matrix') 
    end

end