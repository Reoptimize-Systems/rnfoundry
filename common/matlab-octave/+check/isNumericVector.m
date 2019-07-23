function ok = isNumericVector (inputval, throw, name, posorneg)
% checks if input is a real numeric matrix (m x n), optionally ALL +ve or -ve only
%
% Syntax
%
%  ok = isNumericVector (inputval, throw)
%  ok = isNumericVector (..., name)
%  ok = isNumericVector (..., posorneg)
%
% Input
%
%  inputval - value to be tested if it is a real numeric vector (for which
%   isvector(inputval) returns true, i.e. it has only 2 dimensions).
%
%  throw - logical flag determining whether an error is thrown by
%   isNumericVector if num fails check
%
%  name - optional string used to customise the error message. The error
%   will be "<name> must be a real numeric vector (2D)". Default is 'input'
%   if not supplied.
%
%  posorneg - optional scalar integer flag indicating if an additional
%   check for the number being positive or negative  should be performed.
%   If posorneg == 0, no additional check is performed. If posorneg > 0,
%   isNumericVector checks if all values in 'num' are greater than or equal
%   to zero. If posorneg < 0, isNumericVector checks if all values in 'num'
%   are less than or equal to zero. If checks fail ok is set to false, and
%   an error is thrown depending on the value of 'throw'. Default is zero.
%
% Output
%
%  ok - logical flag indicating if check was passed
%
% See also: check.isNumeric, check.isNumericScalar
%

    if nargin < 3
        name = 'input';
    end
    
    if nargin < 4
        posorneg = 0;
    end
    
    ok = true;
    
    if ~( isnumeric (inputval) && isvector (inputval) && isreal (inputval) )
        ok = false;
        if throw
            error ('%s must be a real numeric vector (1D)', name);
        end
    elseif posorneg > 0 && any(inputval < 0)
        ok = false;
        if throw
            error ('%s must be a real numeric vector (1D) with all values greater than or equal to zero');
        end
    elseif posorneg < 0 && any(inputval > 0)
        ok = false;
        if throw
            error ('%s must be a real numeric vector (1D) with all values less than or equal to zero');
        end
    end

end

