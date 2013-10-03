function MaxSigma3 = Table32r1bMaxSigma3(vars)
% Table32r1fMaxSigma2: Calculates the maximum circumferential normal stress
% (sigma_2) in a cylinder undergoing a uniform internal radial pressure
% with the ends capped, as calculated in 'Roark's Formulas for Stress &
% Strain 6th edition' in table 32, page 638 case 1 manner 1b. 
%
% Input: 
%   
%   vars - (n x 1) matrix:-
%          Col 1. q, the unit pressure on the vessel (force per unit area)
%
% Output:
%
%   MaxSigma2 - (n x 1) column vector of values of MaxSigma3, the maximum
%               radial normal stress in the cylinders
%

    if size(vars,2) == 1
        
        q = vars(:,1);
        
        % max sigma_2 occurs when r = b
        MaxSigma3 = -q;
       
    else
       error('Matrix dimensions do not agree. Table32r1bMaxSigma3 requires a (n x 1) column matrix') 
    end


end