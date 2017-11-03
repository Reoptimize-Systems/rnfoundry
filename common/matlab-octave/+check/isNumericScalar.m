function ok = isNumericScalar (num, throw, name)
    % checks if input is a real scalar value
    %
    % Syntax
    %
    %  ok = isNumericScalar (vec, throw)
    %  ok = isNumericScalar (..., name)
    %
    % Input
    %
    %  num - value to be tested if it is a real numeric
    %    scalar
    %
    %  throw - logical flag determining whether an error is thrown
    %   by isNumericScalar if num fails check
    %
    %  name - optional string used to customise the error message.
    %   The error will be <name> must be a scalar numeric value.
    %   Default is 'value' if not supplied.
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
    if ~( isnumeric (num) && isscalar (num) && isreal (num) )

        ok = false;

        if throw
            error ('%s must be a scalar numeric value', name);
        end
    end

end

