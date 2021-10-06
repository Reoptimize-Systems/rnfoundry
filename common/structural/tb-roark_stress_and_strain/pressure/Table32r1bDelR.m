function delR = Table32r1bDelR(vars)
% function: Table32r1bDelR
% 
% Calculates radial displacement of the circumference of a cylinder
% undergoing a uniform internal pressure q in all directions, as calculated
% in 'Roark's Formulas Stress & Strain 6th edition' in table 32, page 638
% case 1 manner 1b. 
%
% Input: 
%   
%   vars - (n x 5) matrix:-
%          Col 1. q, the unit pressure on the vessel (force per unit area)
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
        
        q = vars(:,1);
        a = vars(:,2);
        b = vars(:,3);
        E = vars(:,4);
        v = vars(:,5);
        
        % first we calculate the displacement of the outer radius (delta a)
        delR(:,1) = (q .* a ./ E) .* b.^2 .* (2 - v) ./ (a.^2 - b.^2);
        
        delR(:,2) =  (q .* b ./ E) .* (a.^2 .* (1 + v) + b.^2 .* (1 - 2 .* v)) ./ (a.^2 - b.^2);
        
    else
       error('Matrix dimensions do not agree. Table32r1bDelR requires a (n x 5) column matrix') 
    end

end