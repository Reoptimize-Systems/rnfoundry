function sigma = Table32r1fNormStresses(vars)
% function: Table32r1fNormStresses
% 
% Calculates normal stresses in the longitudinal, circumferential and
% radial directions for a thick walled cylinder with a linearly varying
% radial body force varying from some peak value at the inner radius to
% zero at the outer radius. This function is based on the formulas in Table
% 32, case 1f in "Roark Formulas for Stress and Strain, 6th Edition" page
% 639.
%
% Input: 
%   
%   vars - (n x 5) matrix:-
%          Col 1. ldeltab, the unit pressure on the vessel (force per unit volume)
%          Col 2. a, the outer radius
%          Col 3. b, the inner radius
%          Col 4. v, Poisson's ratio for the material
%          Col 5. r, the radial distnce from the centre at which the
%                 stresses are to be calculated 
%
% Output:
%
%   sigma - (n x 3) column vector of values of sigma, the normal stresses in
%           the longitudinal, circumferential and radial directions
%           respectively.
%           
%

    if size(vars,2) == 5
        
        ldeltab = vars(:,1);
        a = vars(:,2);
        b = vars(:,3);
        v = vars(:,4);
        r = vars(:,5);
        
        % There are no longitudinal normal stresses, so sigma_1 is zero
        sigma(1:size(a,1),1) = 0;
        
        % We will calculate sigma_2, the circumferential normal stresses in
        % stages for clarity
        sigma(:,2) = ((7 + 5.*v).*a.^4 - 8.*(2 + v).*a.*b.^3 + 3.*(3 + v).*b.^4) ./ (24.*(a - b).*(a.^2 - b .^2));
        
        % The first term in sigma_3 is identical to that in sigma_2, so to
        % avoid repitition of calculation we simply copy this first result.
        sigma(:,3) = sigma(:,2);
        
        
        % continue calculating sigma_2
        sigma(:,2) = sigma(:,2) - (1 + 2.*v).*a.*r ./ (3.*(a - b));
        
        sigma(:,2) = sigma(:,2) + (1 + 3.*v).*r.^2 ./ (8.*(a-b));
        
        sigma(:,2) = sigma(:,2) + ((b.^2.*a.^2) ./ (24 .* r.^2)) .* (((7 + 5.*v).*a.^2 - 8.*(2 + v).*a.*b + 3.*(3 + v).*b.^2) ./ ((a - b).*(a.^2 - b .^2)));
        
        sigma(:,2) = sigma(:,2) .* ldeltab;
        
        
        % Continue calculating sigma_3
        sigma(:,3) = sigma(:,3) - ((2 + v).*a.*r) ./ (3.*(a - b));
        
        sigma(:,3) = sigma(:,3) + ((3 + v).* r.^2) ./ ( 8.*(a - b));
        
        sigma(:,3) = sigma(:,3) - ((b.^2.*a.^2) ./ (24 .* r.^2)) .* (((7 + 5.*v).*a - 3.*(3 + v).*b) ./ (a.^2 - b .^2));
        
        sigma(:,3) = sigma(:,3) .* ldeltab;
        
    else
        
       error('Matrix dimensions do not agree. Table32r1fNormStresses requires an (n x 5) matrix') 
       
    end

end