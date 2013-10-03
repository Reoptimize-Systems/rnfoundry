function delR = Table28r1bDelR(vars)
% function: Table28r1bDelR
% 
% Calculates radial displacement of the circumference of a cylinder
% undergoing a unform axial loading, as calculated in 'Roark's Formulas 
% Stress & Strain 6th edition' in table 28, page 518 case 1 manner 1b. 
%
% Input: 
%   
%   vars - (n x 5) matrix:-
%          Col 1. q, the unit pressure on the vessel (force per unit area)
%          Col 2. R, the outer radius
%          Col 3. E, youngs modulus for the material
%          Col 4. t, the annulus (wall) thickness
%
% Output:
%
%   delR - (n x 1) column vector of values of delR, the radial displacement
%   of the circumference of the cylinder
%
    if size(vars,2) == 4
        delR = vars(:,1) .* (vars(:,2).^2) ./ (vars(:,3) .* vars(:,4));
    else
       error('Matrix dimensions do not agree. Table28r1a requires a (n x 4) column matrix') 
    end

end