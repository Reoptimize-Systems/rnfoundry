function y_c = CircleY (vars)
% Calculates the maximum distance from the neutral axis to the extremities
% in the Y direction of a solid circle, as calculated in 'Roark's Formulas
% Stress & Strain'
%
% Input: 
%   
%   IVars - (n x 1) matrix, containing values of R, the radius
%
% Output:
%
%   y_c - (n x 1) column vector of values of y_c, the distance to the
%     extremeties
%

    y_c = vars (:,1);
    
end