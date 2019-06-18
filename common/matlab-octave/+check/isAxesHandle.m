function ok = isAxesHandle (input, throw, inputname)
% checks if input is an axes handle
%
% Syntax
%
%  ok = isAxesHandle (input, throw)
%  ok = isAxesHandle (..., inputname)
%
% Input
%
%  input - value to be tested if it is of type axes handle
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
    
    if isoctave ()
        % do nothing
        ok = isaxes (input);
    else
        ok = isa (input, 'matlab.graphics.axis.Axes');
    end
    
    if throw && ~ok

        error ('%s must be an axes handle (or array of axes handles)', inputname);

    end

end