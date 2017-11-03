function ok = isLogicalScalar (tfvalue, throw, name)
    % checks if input is a scalar logical value
    %
    % Syntax
    %
    %  ok = isLogicalScalar (tfvalue, throw)
    %  ok = isLogicalScalar (..., name)
    %
    % Input
    %
    %  tfvalue - value to be tested if it is a scalar logical value
    %
    %  throw - logical flag determining whether an error is thrown
    %   by isLogicalScalar if tfvalue fails check
    %
    %  name - optional string used to customise the error message.
    %   The error will be "<name> must be a scalar logical
    %   (true/false) value". Default is 'value' if not supplied.
    %
    % Output
    %
    %  ok - logical flag indicating if check was passed
    %

    if nargin < 3
        name = 'value';
    end

    ok = true;
    if ~( islogical (tfvalue) && isscalar (tfvalue) )

        ok = false;

        if throw
            error ('%s must be a scalar logical value', name);
        end
    end

end