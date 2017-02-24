classdef const < mbdyn.pre.driveCaller
    
    properties (GetAccess = public, SetAccess = private)
        value;
    end
    
    methods
        
        function self = const (value)
            
            assert ( isscalar (value) && isnumeric (value), ...
                'value must be a numeric scalar');
            
            self.value = value;
            
        end
        
        function str = generateOutputString (self)
            
            str = sprintf ('const, %g', self.value);
            
        end
        
    end
    
end