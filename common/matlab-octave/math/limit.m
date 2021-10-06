function M = limit (M, range)
% limit a value or values in a matrix to a given range
%
% Syntax
%
% limitedM = limit (M, range)
%
% Input
%
%   M - matrix of values to be restricted in size
%
%   range - two element vector of values, the first value is the lower
%     limit to be placed on the values in M, the second is the upper limit
%     to be placed on the values in M
%
% Output
%
%   limitedM - Input matrix M, but with all values less than range(1)
%     replaced with range(1) and all values greater than range(2) replaced
%     with range(2)
%

% Created by Richard Crozier 2014


    if numel (range) ~= 2  || ~isnumeric (range)
        error ('range must be a 1 or 2 element numeric vector.') 
    end
    
    if range(1) >= range(2)
        error ('range(1) must be less than range(2).')
    end

    M(M < range(1)) = range(1);
    
    M(M > range(2)) = range(2);
    
end