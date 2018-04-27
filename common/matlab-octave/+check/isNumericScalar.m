function ok = isNumericScalar (num, throw, name, posorneg)
% checks if input is a real scalar value, optionally +ve or -ve only
%
% Syntax
%
%  ok = isNumericScalar (vec, throw)
%  ok = isNumericScalar (..., name)
%  ok = isNumericScalar (..., posorneg)
%
% Input
%
%  num - input to be tested if it is a real numeric scalar
%
%  throw - logical flag determining whether an error is thrown by
%   isNumericScalar if num fails check
%
%  name - optional string used to customise the error message. The error
%   will be <name> must be a scalar numeric value. Default is 'value' if
%   not supplied.
%
%  posorneg - optional scalar integer flag indicating if an additional
%   check for the number being positive or negative  should be performed.
%   If posorneg == 0, no additional check is performed. If posorneg > 0,
%   isNumericScalar checks if 'num' is greater than or equal to zero. If
%   posorneg < 0, isNumericScalar checks if 'num' is less than or equal to
%   zero. If checks fail ok is set to false, and an error is thrown
%   depending on the value of 'throw'. Default is zero.
%
% Output
%
%  ok - logical flag indicating if check was passed
%
%
% See also: check.isNumeric, check.isNumericMatrix
%

    if nargin < 3
        name = 'value';
    end
    
    if nargin < 4
        posorneg = 0;
    end
    
    ok = true;
    
    if ~( isnumeric (num) && isscalar (num) && isreal (num) )
        ok = false;
        if throw
            error ('%s must be a real scalar numeric value', name);
        end
    elseif posorneg > 0 && num < 0
        ok = false;
        if throw
            error ('%s must be a real scalar numeric value greater than or equal to zero');
        end
    elseif posorneg < 0 && num > 0
        ok = false;
        if throw
            error ('%s must be a real scalar numeric value less than or equal to zero');
        end
    end

end

