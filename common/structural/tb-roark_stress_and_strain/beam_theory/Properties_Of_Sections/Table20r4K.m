function K = Table20r4K(vars)
% Table20r4K: calculates the torsion constant for a solid rectangular
% cross-section according to the formulas in Table 20 row 4 (and 3 for
% square cross-sections) of Roark's Formulas for Stress and Strain, 6th
% Edition
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
%   IVars - (n x 2) matrix of values as below
%           IVars(:,1) - d,rectangle height
%           IVars(:,2) - b, width
%
% Output:
%
%   K - (n x 1) column vector of values of K, the torsion constant of a
%   beam with a solid rectangular cross-section
%

    d = vars(:,1);
    b = vars(:,2);

    K(b > d, 1) = ((b(b>=d)./2)) .* ((d(b>=d)./2))^3 .* ...
        ((16/3) - 3.36 .* ((d(b>=d)./2)./(b(b>=d)./2)) .* (1 - ((d(b>=d)./2).^4 ./ (12.*(b(b>=d)./2).^4))));
    
    K(b < d, 1) = ((d(b<d)./2)) .* ((b(b<d)./2))^3 .* ...
        ((16/3) - 3.36 .* ((b(b<d)./2) ./ (d(b<d)./2)) .* (1 - ((b(b<d)./2).^4 ./ (12.*(d(b<d)./2).^4))));
    
    K(b == d) = 2.25 .* (b(b==d)./2).^4;

end