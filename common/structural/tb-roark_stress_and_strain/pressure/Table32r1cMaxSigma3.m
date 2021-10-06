function MaxSigma3 = Table32r1cMaxSigma3(vars)
% Table32r1fMaxSigma3: Calculates the maximum radial normal stress
% (sigma_3) in a cylinder undergoing a uniform external radial pressure
% with longitudinal pressure zero or balanced externally, as calculated in
% 'Roark's Formulas for Stress & Strain 6th edition' in table 32, page 638
% case 1 manner 1c. 
%
% Input: 
%   
%   vars - (n x 1) matrix:-
%          Col 1. q, the unit pressure on the vessel (force per unit area)
%
% Output:
%
%   MaxSigma2 - (n x 1) column vector of values of MaxSigma3, the maximum
%               shear stress in the cylinders
%

    if size(vars,2) == 1
        
        q = vars(:,1);
        
        % max sigma_2 occurs when r = b
        MaxSigma3 = -q;
       
    else
       error('Matrix dimensions do not agree. Table32r1cMaxShearStress requires a (n x 1) column matrix') 
    end


end