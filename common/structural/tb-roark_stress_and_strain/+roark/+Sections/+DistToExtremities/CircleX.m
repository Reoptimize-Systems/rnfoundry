function x_c = CircleX (vars)
% Calculates the maximum distance from the neutral axis to the extremities
% in the X direction of a cirlce, as calculated in 'Roark's Formulas Stress
% & Strain'
%
% Input: 
%   
%   IVars - (n x 1) matrix, containing values of R, the
%     radius
%
% Output:
%
%   x_c - (n x 1) column vector of values of x_c, the distance to the
%     extremeties
%

    x_c = vars (:,1);
    
end