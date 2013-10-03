function MaxSigma2 = Table32r1eMaxSigma2(vars)
% Table32r1eMaxShearStress: Calculates the maximum longitudinal stresses in a
% cylinder undergoing a uniform radial body force, as calculated in
% 'Roark's Formulas for Stress & Strain 6th edition' in table 32, page 638
% case 1 manner 1e. 
%
% Input: 
%   
%   vars - (n x 4) matrix:-
%          Col 1. delta, the unit pressure on the vessel (force per unit area)
%          Col 2. a, the outer radius
%          Col 3. b, the inner radius
%          Col 4. v, Poisson's ratio for the material
%
% Output:
%
%   MaxSigma - (n x 1) column vector of values of MaxSigma2, the maimum
%              circumferential normal stress in the cylinders
%

    if size(vars,2) == 4
        
        delta = vars(:,1);
        a = vars(:,2);
        b = vars(:,3);
        v = vars(:,4);
        
        % max sigma_2 occurs at r = b
        MaxSigma2 = (delta .* a.^2 ./ 3) .* (2 .* (2 + v) ./ (a + b) + b .* (1 - v) ./ a.^2);
       
    else
       error('Matrix dimensions do not agree. Table32r1eMaxShearStress requires a (n x 4) column matrix') 
    end
end