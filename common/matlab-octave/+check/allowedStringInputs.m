function [ ok, allowed_ind ] = allowedStringInputs (input, allowedstrs, throw, inputname)
% checks if string is one of a set of allowed values
%
% Syntax
%
%  ok = allowedStringInputs (input, allowedstrs, throw)
%  ok = allowedStringInputs (..., inputname)
%
% Input
%
%  input - value to be tested if it is a valid string. The string will be
%    compared to the list of valid strings in allowedstrs.
%
%  allowedstrs - cell string array of valid strings for comparison with
%    the string provided in input.
%
%  throw - logical flag determining whether an error is thrown by
%    checkOrientationDescription if input fails check
%
%  inputname - (optional) string with name to use in error thrown by
%    checkAllowedStringInputs (if 'throw' is true). The error will be
%    "<name> must be one of: 'string 1' | 'string 2' | 'string 3'" where
%    allowedstrs is {'string 1', 'string 2', 'string 3'}.
%
% Output
%
%  ok - logical flag indicating if check was passed
%

    if nargin < 4
        inputname = 'input';
    end

    if ~iscellstr (allowedstrs)
        error ('allowedstrs must be a cell array of strings containing a list of allowed values for the input');
    end

    ok = true;
    
    allowed_ind = find (strcmp (input, allowedstrs));

    if isempty(allowed_ind)
        ok = false;
    end

    if throw && ~ok
        
        if numel (allowedstrs) == 1
            error ('%s must be ''%s'' ', inputname, allowedstrs{1});
        else
            error ('%s must be one of: %s ', inputname, sprintf (['''%s''', repmat(' | ''%s''', 1, numel(allowedstrs)-1)], allowedstrs{:}));
        end

    end

end
