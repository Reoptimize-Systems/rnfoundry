function I = Table1r2I_1(IVars)
% Table1r2I_1: Calculates the second moment of inertia of a beam 
% with a rectangular cross-section, as calculated in 'Roark's 
% Formulas Stress & Strain 6th edition' in table 1, page 62 row 2.
%
%             2
%             .
%        _____._______
%    ^  |     .       |
%    |  |     .       |
%    |  |     .       |
% d  |  |.............|......... 1
%    |  |     .       |
%    |  |     .       |
%    |  |     .       |
%    v  |_____._______|
%       <-----.------>
%          b  .
%             .
%
% Input: 
%
%   IVars - (n x 4) matrix of values as below
%           IVars(:,1) - d,rectangle height
%           IVars(:,2) - b, width
%
% Output:
%
%   I - (n x 1) column vector of values of I, the moment of inertia for a
%       beam with a rectangular cross-section
% 
    if size(IVars,2) == 2
        d = IVars(:,1);
        b = IVars(:,2);

        I = (d.^3) .* b ./ 12;
    else
        error('Dimensions of Input variables incorrect. Table1r2I_2 accepts a (n x 2) column vector of values')
    end
     
end