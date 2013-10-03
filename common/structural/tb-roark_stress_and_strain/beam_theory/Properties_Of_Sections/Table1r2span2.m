function span2val = Table1r2span2(vars)
% Table1r2span1: Calculates the maximum span of a beam with a
% rectangular cross-section, in the axis 1
% 
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
%   outer radius, the second contains values of t, the annulus thickness.
%
% Output:
%
%   span2val - (n x 1) column vector of values of the maximum span of the
%   section in axis 2

    span2val = vars(:,1); 
    
end