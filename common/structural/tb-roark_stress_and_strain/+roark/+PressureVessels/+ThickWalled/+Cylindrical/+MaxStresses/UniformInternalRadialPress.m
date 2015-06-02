function MaxShearStress = UniformInternalRadialPress (vars)
% Calculates the maximum shear stress for a thick walled cylinder with a
% uniform internal radial pressure for a disk or a shell based on the
% formulas in "Roark's Formulas for Stress and Strain".
%
% Input: 
%   
%   vars - (n x 3) matrix:-
%    Col 1. q, the unit pressure on the vessel (force per unit area)
%    Col 2. a, the outer radius
%    Col 3. b, the inner radius
%
% Output:
%
%   MaxSigma2 - (n x 1) column vector of values of the maximum shear stress
%     in the cylinders
%

    if size(vars,2) == 3
        
        q = vars(:,1);
        a = vars(:,2);
        b = vars(:,3);
        
        MaxShearStress =  q .* a.^2 ./ (a.^2 - b.^2);
       
    else
       error('Matrix dimensions do not agree. UniformInternalRadialPress requires a (n x 3) column matrix') 
    end
end