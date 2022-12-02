function ok = isScalarInteger (num, throw, name)
% checks if input is a scalar integer value to machine precision
%
% Syntax
%
%  ok = isScalarInteger (num, throw)
%  ok = isScalarInteger (..., name)
%
% Input
%
%  num - value to be tested if it is a scalar integer
%
%  throw - logical flag determining whether an error is thrown
%   by isScalarInteger if num fails check
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
    if ~( isnumeric (num) && isscalar (num) && isreal (num) && check.isint2eps (num) )

        ok = false;

        if throw
            error ('%s must be a scalar integer (to machine precision)', name);
        end
    end

end

