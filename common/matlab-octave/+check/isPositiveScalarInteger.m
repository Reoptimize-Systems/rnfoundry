function ok = isPositiveScalarInteger (num, throw, name)
    % checks if input is a positive scalar integer to machine precision
    %
    % Syntax
    %
    %  ok = isPositiveScalarInteger (num, throw)
    %  ok = isPositiveScalarInteger (..., name)
    %
    % Input
    %
    %  num - value to be tested if it is a scalar integer
    %
    %  throw - logical flag determining whether an error is thrown
    %   by isPositiveScalarInteger if num fails check
    %
    %  name - optional string used to customise the error message.
    %   The error will be "<name> must be a positive ( > 0 ) scalar integer
    %   (to machine precision)". Default is 'value' if not supplied.
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
    if ~( isnumeric (num) ...
            && isscalar (num) ...
            && isreal (num) ...
            && check.isint2eps (num) ...
            && num > 0 )

        ok = false;

        if throw
            error ('%s must be a positive ( > 0 ) scalar integer (to machine precision)', name);
        end
    end

end

