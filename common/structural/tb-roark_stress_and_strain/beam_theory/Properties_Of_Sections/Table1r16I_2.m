function I = Table1r16I_2(IVars)
% function: Table1r16I_2
% 
% Calculates the moment of inertial of a beam with a hollow circular
% cross-section, as calculated in 'Roark's Formulas Stress & Strain 6th
% edition' in table 1, page 66 row 16. This funtion is suited for
% hollow circles with a very thin annulus.
%
% Input: 
%   
%   IVars - (n x 2) matrix, the first column contains values of R, the
%   outer radius, the second contains values of t, the annulus thickness.
%
% Output:
%
%   I - (n x 1) column vector of values of I, the moment of inertia for a
%   beam with hollow circular cross-section and very thin annulus
%
    
    % Same as I1 for circle
    I = Table1r16I_1(IVars);
     
end