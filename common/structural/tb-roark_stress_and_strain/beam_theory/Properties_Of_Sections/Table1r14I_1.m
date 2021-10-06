function I = Table1r14I_1(IVars)
% Table1r14I_1: Calculates the moment of inertia of a beam with a solid
% circular cross-section, as calculated in 'Roark's Formulas Stress &
% Strain 6th edition' in table 1, page 66 row 14.
%
% Input: 
%   
%   IVars - (n x 1) column vector of values of R, the radius of the
%   circular cross-section
%
% Output:
%
%   I - (n x 1) column vector of values of I, the moment of inertia for a
%   beam with solid circular cross-section
%
    if size(IVars,2) > 1
        error('IVars has too many columns, IVars must be a (n x 1) column vector')
    end
    
    if size(IVars,2) == 1
        I = (pi / 4) .* (Vars.^4);
    else
        error('Dimensions of Input variables incorrect. Table1r14I_1 accepts a (n x 1) column vector of values of R')
    end
     
end