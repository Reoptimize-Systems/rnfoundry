classdef joint < mbdyn.pre.element
    
    properties
        
    end
    
    methods
        
        function str = generateMBDynInputString (self)
            str = sprintf ('    joint : %d, %s,', self.label, self.type);
        end
        
    end
    
    methods (Access = protected)
        

        
    end
    
end