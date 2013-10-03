function I = Table1r16I_1(IVars)
% function: Table1r16I_1
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
    if size(IVars,2) > 2
        error('IVars has too many columns, IVars must be a (n x 2) matrix')
    end
    
    if size(IVars,2) == 2
        I = pi .* (IVars(:,1).^3) .* IVars(:,2);
    else
        error('Dimensions of Input variables incorrect. Table1r16I_1 accepts a (n x 2) column vector of values of R and t')
    end
     
end