classdef const < mbdyn.pre.drive
    
    properties (GetAccess = public, SetAccess = private)
        value;
    end
    
    methods
        
        function self = const (value)
            % const drive constructor
            %
            % Syntax
            %
            % cd = mbdyn.pre.const (value)
            %
            % Description
            %
            % Drive which provides a constant value.
            %
            % Input
            %
            %  value - scalar numeric value of the constant
            %
            % Output
            %
            %  cd - mbdyn.pre.const object
            %
            %
            %
            % See Also: 
            %
            
            assert ( isscalar (value) && isnumeric (value), ...
                'value must be a numeric scalar');
            
            self.value = value;
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = sprintf ('const, %s', self.formatNumber (self.value));
            
        end
        
    end
    
end