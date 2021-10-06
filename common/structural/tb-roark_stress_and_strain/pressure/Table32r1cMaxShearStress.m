function MaxShearStress = Table32r1cMaxShearStress(vars)
% Table32r1fMaxShearStress: Calculates the maximum shear stress in a
% cylinder undergoing a uniform external radial pressure with longitudinal
% pressure zero or balanced externally, as calculated in 'Roark's Formulas
% for Stress & Strain 6th edition' in table 32, page 638 case 1 manner 1c. 
%
% Input: 
%   
%   vars - (n x 3) matrix:-
%          Col 1. q, the unit pressure on the vessel (force per unit area)
%          Col 2. a, the outer radius
%          Col 3. b, the inner radius
%
% Output:
%
%   MaxSigma2 - (n x 1) column vector of values of the maximum
%               shear stress in the cylinders
%

    if size(vars,2) == 3
        
        q = vars(:,1);
        a = vars(:,2);
        b = vars(:,3);
        
        MaxShearStress =  q .* a.^2 ./ (a.^2 - b.^2);
       
    else
       error('Matrix dimensions do not agree. Table32r1cMaxShearStress requires a (n x 3) column matrix') 
    end
end