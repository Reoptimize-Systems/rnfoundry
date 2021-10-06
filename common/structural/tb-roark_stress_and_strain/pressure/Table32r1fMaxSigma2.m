function MaxSigma2 = Table32r1fMaxSigma2(vars)
% Table32r1fMaxSigma2: Calculates the maximum circumferential normal
% stress in a cylinder undergoing a linearly varying radial body force, as
% calculated in 'Roark's Formulas for Stress & Strain 6th edition' in table
% 32, page 639 case 1 manner 1f. 
%
% Input: 
%   
%   vars - (n x 4) matrix:-
%          Col 1. deltab, the unit pressure on the vessel (force per unit area)
%          Col 2. a, the outer radius
%          Col 3. b, the inner radius
%          Col 4. v, Poisson's ratio for the material
%
% Output:
%
%   MaxSigma2 - (n x 1) column vector of values of MaxSigma2, the maimum
%               shear stress in the cylinders
%
    
    if size(vars,2) == 4
        
        deltab = vars(:,1);
        a = vars(:,2);
        b = vars(:,3);
        v = vars(:,4);
        
        % max sigma_2 occurs when r = b
        MaxSigma2 = (deltab./12).*((2.*a.^4 + (1+v).* a.^2 .* (5.*a.^2 - 12.*a.*b + 6.*b.^2) -  (1-v).*b.^3 .* (4.*a - 3.*b))./((a - b).*(a.^2 - b.^2)));
       
    else
       error('Matrix dimensions do not agree. Table32r1eMaxShearStress requires a (n x 4) column matrix') 
    end

end