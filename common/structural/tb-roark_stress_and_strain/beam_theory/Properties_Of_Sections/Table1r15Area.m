function I = Table1r15Area(IVars)
% Table1r15Area: Calculates the cross-sectional area of a beam with a
% hollow circular cross-section, as calculated in 'Roark's Formulas Stress
% & Strain 6th edition' in table 1, page 66 row 15. This funtion is not
% suited for hollow circles with a very thin annulus.
%
% Input: 
%   
%   IVars - (n x 2) matrix, the first column contains values of R, the
%   outer radius, the second contains values of Ri, the inner radius.
%
% Output:
%
%   I - (n x 1) column vector of values of I, the moment of inertia for a
%   beam with hollow circular cross-section
%
    if size(IVars,2) > 2
        error('IVars has too many columns, IVars must be a (n x 2) matrix')
    end
    
    if size(IVars,2) == 2
        
        R = IVars(:,1);
        
        Ri = IVars(:,2);
        
        I = (pi / 4) .* (R.^4 - Ri.^4);
    else
        error('Dimensions of Input variables incorrect. Table1r14I_1 accepts a (n x 2) column vector of values of R')
    end
     
end