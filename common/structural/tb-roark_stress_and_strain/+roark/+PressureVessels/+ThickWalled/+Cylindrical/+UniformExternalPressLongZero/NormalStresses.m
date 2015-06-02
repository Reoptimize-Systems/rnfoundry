function sigma = UniformExternalPressLongZero (vars)
% Calculates normal stresses in the longitudinal, circumferential and
% radial directions for a thick walled cylinder with a a uniform external
% radial pressure with longitudinal pressure zero or balanced externally.
% This function is based on the formulas in "Roark's Formulas for Stress
% and Strain".
%
% Input: 
%   
%   vars - (n x 4) matrix:-
%          Col 1. q, the unit pressure on the vessel (force per unit area)
%          Col 2. a, the outer radius
%          Col 3. b, the inner radius
%          Col 4. r, the radial distnce from the centre at which the
%                 stresses are to be calculated 
%
% Output:
%
%   sigma - (n x 3) column vector of values of sigma, the normal stresses in
%           the longitudinal, circumferential and radial directions
%           respectively.
%           
%

    if size(vars,2) == 4
        
        q = vars(:,1);
        a = vars(:,2);
        b = vars(:,3);
        r = vars(:,4);
        
        % There are no longitudinal normal stresses, so sigma_1 is zero
        sigma(1:size(a,1),1) = 0;
        
        % Sigma_2, the circumferential normal stresses
        sigma(:,2) = -q .* a.^2 .* (b.^2 + r.^2) ./ (r.^2 .* (a.^2 - b.^2));
        
        % Sigma_3
        sigma(:,3) = -q .* a.^2 .* (r.^2 - b.^2) ./ (r.^2 .* (a.^2 - b.^2));
        
    else
        
       error('Matrix dimensions do not agree. Table32r1fNormStresses requires an (n x 4) matrix') 
       
    end

end