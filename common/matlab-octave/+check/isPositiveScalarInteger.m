function ok = isPositiveScalarInteger (num, throw, name, allowzero)
% checks if input is a real positive scalar integer to machine precision
%
% Syntax
%
% ok = isPositiveScalarInteger (num, throw)
% ok = isPositiveScalarInteger (..., name)
% ok = isPositiveScalarInteger (..., allowzero)
%
% Input
%
%  num - value to be tested if it is a positive scalar integer. By default
%   it must also be greater than zero, this behaviour can be changed with
%   the allowzero option.
%
%  throw - logical flag determining whether an error is thrown by
%   isPositiveScalarInteger if num fails check
%
%  name - optional string used to customise the error message. The error
%   will be "<name> must be a positive ( > 0 ) scalar integer (to machine
%   precision)" or "<name> must be a positive ( >= 0 ) scalar integer (to
%   machine precision)". Default is 'value' if not supplied.
%
%  allowzero - opional true/false flag indicating whether to allow
%   integers >= 0 rather than only integers > 0 (the default).
%
% Output
%
%  ok - logical flag indicating if check was passed
%

    if nargin < 3
        name = 'value';
    end
    
    if nargin < 4
        allowzero = false;
    end
    
    check.isLogicalScalar (allowzero, true, 'allowzero');

    ok = true;
    
    if allowzero
        
        if ~( isnumeric (num) ...
            && isscalar (num) ...
            && isreal (num) ...
            && check.isint2eps (num) ...
            && num >= 0 )
    
            ok = false;

            if throw
                error ('%s must be a positive ( > 0 ) scalar integer (to machine precision)', name);
            end
        
        end
        
    else
        
        if ~( isnumeric (num) ...
            && isscalar (num) ...
            && isreal (num) ...
            && check.isint2eps (num) ...
            && num > 0 )
    
            ok = false;

            if throw
                error ('%s must be a positive ( >= 0 ) scalar integer (to machine precision)', name);
            end
        
        end
        
    end

end

