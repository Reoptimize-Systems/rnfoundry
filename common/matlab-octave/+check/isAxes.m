function ok = isAxes (obj, throw, name)
% checks if input is an axes object (or axes handle)
%
% Syntax
%
%  ok = isAxes (obj, throw)
%  ok = isAxes (..., name)
%
% Input
%
%  obj - value to be tested if it is an axes object or handle
%
%  throw - logical flag determining whether an error is thrown
%   by checkIsAxes if obj fails check
%
%  name - optional string used to customise the error message.
%   The error will be <name> must be a scalar integer (to
%   machine precision). Default is 'value' if not supplied.
%
%
% Output
%
%  ok - logical flag indicating if check was passed
%

    if nargin < 3
        name = 'value';
    end

    ok = true;
    if isoctave ()
        if ~all (isaxes (obj))
            ok = false;
        end
    else
        if ~isa (obj, 'matlab.graphics.axis.Axes')
            ok = false;
        end
    end

    if throw && ~ok
        error ('%s must be a scalar integer (to machine precision)', name);
    end

end