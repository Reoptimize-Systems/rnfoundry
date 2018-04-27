function ok = isNumeric (inputval, throw, name, posorneg)
% checks if input is a real and numeric, optionally ALL +ve or -ve only
%
% Syntax
%
%  ok = isNumeric (inputval, throw)
%  ok = isNumeric (..., name)
%  ok = isNumeric (..., posorneg)
%
% Input
%
%  inputval - input to be tested if it is composed of real numeric values
%   (can be a matrix of any size and number of dimensions).
%
%  throw - logical flag determining whether an error is thrown by
%   isNumeric if num fails check
%
%  name - optional string used to customise the error message. The error
%   will be "<name> must be all real numeric values". Default is 'input' if
%   not supplied.
%
%  posorneg - optional scalar integer flag indicating if an additional
%   check for the number being positive or negative  should be performed.
%   If posorneg == 0, no additional check is performed. If posorneg > 0,
%   isNumeric checks if all values in 'num' are greater than or equal to
%   zero. If posorneg < 0, isNumeric checks if all values in 'num' are less
%   than or equal to zero. If checks fail ok is set to false, and an error
%   is thrown depending on the value of 'throw'. Default is zero.
%
% Output
%
%  ok - logical flag indicating if check was passed
%
%
% See also: check.isNumericMatrix, check.isNumericScalar
%

    if nargin < 3
        name = 'input';
    end
    
    if nargin < 4
        posorneg = 0;
    end
    
    ok = true;
    
    if ~( isnumeric (inputval) && isreal (inputval) )
        ok = false;
        if throw
            error ('%s must be all real numeric values', name);
        end
    elseif posorneg > 0 && any(inputval < 0)
        ok = false;
        if throw
            error ('%s must be all real numeric values greater than or equal to zero');
        end
    elseif posorneg < 0 && any(inputval > 0)
        ok = false;
        if throw
            error ('%s must be a all real numeric values less than or equal to zero');
        end
    end

end

