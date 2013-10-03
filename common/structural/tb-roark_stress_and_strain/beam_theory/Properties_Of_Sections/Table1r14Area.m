function A = Table1r14Area(vars)
% function: Table1r14Area: Calculates the cross-sectional area of a beam
% with a solid circular cross-section, as calculated in 'Roark's Formulas
% Stress & Strain 6th edition' in table 1, page 66 row 14.
%
% Input: 
%   
%   IVars - (n x 1) column vector of values of R, the radius of the
%   circular cross-section
%
% Output:
%
%   A - (n x 1) column vector of values of A, the cross-sectional area of a
%   beam with solid circular cross-section
%
    if size(IVars,2) > 1
        error('IVars has too many columns, IVars must be a (n x 1) column vector')
    end
    
    if size(IVars,2) == 1
        
        r = vars(:,1);
        
        A = pi .* (r.^2);
        
    else
        error('Dimensions of Input variables incorrect. Table1r14I_1 accepts a (n x 1) column vector of values of R')
    end
     
end