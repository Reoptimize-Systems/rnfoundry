function x_c = HollowCircleX (vars)
% Calculates the maximum distance from the neutral axis to the extremities
% in the X direction, as calculated in 'Roark's Formulas Stress & Strain'
%
% Input: 
%   
%   IVars - (n x 2) matrix, the first column contains values of R, the
%     outer radius, the second contains values of Ri, the inner radius.
%
% Output:
%
%   x_c - (n x 1) column vector of values of x_c, the distance to the
%     extremeties
%

    x_c = vars (:,1);
    
end