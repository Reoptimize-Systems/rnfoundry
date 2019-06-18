function ok = isCharOrScalarString (input, throw, inputname)
% checks if input is a character vector or a single string
%
% Syntax
%
%  ok = isCharOrScalarString (input, throw)
%  ok = isCharOrScalarString (..., inputname)
%
% Input
%
%  input - value to be tested if it is a character vecotr os a scalar
%    string array
%
%  throw - logical flag determining whether an error is thrown by
%    isCharOrScalarString if input fails check
%
%  inputname - (optional) string with name to use in error thrown by
%    structHasAllFields (if 'throw' is true). The error will be either
%    "<name> must be a character vector or a scalar string array" or
%    "<name> must be a scalar sting, with numel(<name>) == 1" if the input
%    is an array of string objects with more than one element.
%
% Output
%
%  ok - logical flag indicating if check was passed
%

    if nargin < 3
        inputname = 'input';
    end
    
    ok = true;
    
    if ischar (input)
        % do nothing
    elseif isstring (input)
        % check it's scalar
        if numel (input) ~= 1
            ok = false;
            if throw
                error ('%s must be a scalar sting, i.e. with numel(%s) == 1', inputname, inputname);
            end
        end
    else
         ok = false;
    end
    
    if throw && ~ok

        error ('%s must be a character vector or a scalar string array', inputname);

    end

end