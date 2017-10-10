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