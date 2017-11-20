classdef base < handle
    % base class for the nemoh preprocessing tools
    
    properties (GetAccess = public, SetAccess = protected)
        
    end
    
    methods (Static)
        
        function out = makeCellIfNot (in)
            if ~iscell (in)
                out = {in};
            else
                out = in;
            end
        end
        
        
        function ok = checkNumericScalar (num, throw, name)
            % checks if input is a real scalar value
            %
            % Syntax
            %
            %  ok = checkNumericScalar (vec, throw)
            %  ok = checkNumericScalar (..., name)
            %
            % Input
            %
            %  num - value to be tested if it is a real numeric
            %    scalar
            %
            %  throw - logical flag determining whether an error is thrown
            %   by checkNumericScalar if num fails check
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
        
        function ok = checkLogicalScalar (tfvalue, throw, name)
            % checks if input is a scalar logical value
            %
            % Syntax
            %
            %  ok = checkLogicalScalar (tfvalue, throw)
            %  ok = checkLogicalScalar (..., name)
            %
            % Input
            %
            %  tfvalue - value to be tested if it is a scalar logical value
            %
            %  throw - logical flag determining whether an error is thrown
            %   by checkLogicalScalar if tfvalue fails check
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
        
        function ok = checkScalarInteger (num, throw, name)
            % checks if input is a scalar integer value to machine precision
            %
            % Syntax
            %
            %  ok = checkScalarInteger (num, throw)
            %  ok = checkScalarInteger (..., name)
            %
            % Input
            %
            %  num - value to be tested if it is a scalar integer
            %
            %  throw - logical flag determining whether an error is thrown
            %   by checkScalarInteger if num fails check
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
            if ~( isnumeric (num) && isscalar (num) && isreal (num) && isint2eps (num) )
                
                ok = false;
                
                if throw
                    error ('%s must be a scalar integer (to machine precision)', name);
                end
            end
            
        end
        
        function ok = checkIsAxes (obj, throw, name)
            % checks if input is an axes object (or axes handle)
            %
            % Syntax
            %
            %  ok = checkIsAxes (obj, throw)
            %  ok = checkIsAxes (..., name)
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
                if ~isaxes (obj)
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
        
        function numstr = formatNumber (num)
            % fomats a decimal number for pretty output to nemoh file
            %
            % If it is an integer (to machine precision), it will be output
            % with one decimal place. If a float it will be output to 14
            % decimal places, but with any trailing zeros stripped to
            % shorten it.
            
            if isint2eps (num)
                numstr = sprintf ('%d.0', num);
            else
                numstr = sprintf ('%.18f', num);
                % strip trailing zeros from decimals
                n = numel (numstr);
                while numstr(n) ~= '.'
                    if numstr(n) == '0'
                        numstr(n) = [];
                    else
                        break;
                    end
                    n = n - 1;
                end
            end
            
        end
        
        function numstr = formatInteger (num)
            % fomats a non-decimal integer number for output to nemoh file
            %
            
            if isint2eps (num)
                numstr = sprintf ('%d', num);
            else
                error ('Supplied number is not an integer (to machine precision), cannot format');
            end
            
        end
        
    end
    
end



function result = isint2eps(X)
% isint2eps determines if the numbers in a matrix are integers to the limit
% of the machine floating-point relative accuracy for those values

    theMod = mod(X,1);
    
    result = theMod <= eps(theMod);

end


function t = isoctave()
% ISOCTAVE.M
% ISOCTAVE  True if the operating environment is octave.
%    Usage: t=isoctave();
% 
%    Returns 1 if the operating environment is octave, otherwise
%    0 (Matlab)
% 
% ---------------------------------------------------------------
%
% COPYRIGHT : (c) NUHAG, Dept.Math., University of Vienna, AUSTRIA
%             http://nuhag.eu/
%             Permission is granted to modify and re-distribute this
%             code in any manner as long as this notice is preserved.
%             All standard disclaimers apply.

    if exist('OCTAVE_VERSION')
        % Only Octave has this variable.
        t=1;
    else
        t=0;
    end

end