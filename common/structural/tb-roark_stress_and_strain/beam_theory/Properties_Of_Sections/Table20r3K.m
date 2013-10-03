function K = Table20r3K(vars)
% Table20r3K: calculates the torsion constant for a solid square
% cross-section according to the formula in Table 20 row 3 of Roark's
% Formulas for Stress and Strain, 6th Edition
%
%             2
%             .
%        _____._______
%    ^  |     .       |
%    |  |     .       |
% b  |  |.............|......... 1
%    |  |     .       |
%    |  |     .       |
%    v  |_____._______|
%       <-----.------>
%          b  .
%             .
%  
% Input: 
%   
%   IVars - (n x 1) matrix of values as below
%           IVars(:,1) - b, side length
%
% Output:
%
%   K - (n x 1) column vector of values of K, the torsion constant of a
%   beam with a solid square cross-section
%
    b = vars(:,1);

    K = 2.25 .* (b./2).^4;

end