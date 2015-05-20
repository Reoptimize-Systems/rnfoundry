function y_c = HollowCircleY (vars)
% Calculates the maximum distance from the neutral axis to the extremities
% in the Y direction, as calculated in 'Roark's Formulas Stress & Strain'
%
% Input: 
%   
%   IVars - (n x 2) matrix, the first column contains values of R, the
%   outer radius, the second contains values of Ri, the inner radius.
%
% Output:
%
%   y_c - (n x 1) column vector of values of y_c, the distance to the
%     extremeties
%

    y_c = vars (:,1);
    
end