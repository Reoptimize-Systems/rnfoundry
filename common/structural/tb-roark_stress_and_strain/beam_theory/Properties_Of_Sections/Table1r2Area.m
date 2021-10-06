function A = Table1r2Area(vars)
% Table1r2Area: Calculates the cross-sectional area of a beam with a
% rectangular cross-section, (wide-flange beam with equal flanges) as
% calculated in 'Roark's Formulas Stress & Strain 6th edition' in table 1,
% page 64 row 2. 
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
%   A - (n x 1) column vector of values of A, the cross-sectional area of a
%   beam with a rectangular cross-section
%

    d = vars(:,1);
    b = vars(:,2);
    
    A = b .* d;

end