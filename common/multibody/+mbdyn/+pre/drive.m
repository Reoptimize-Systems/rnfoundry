classdef drive < mbdyn.pre.base
    
    properties
        
    end
    
    methods
        
        function self = drive ()
            
            
        end
        
        function str = generateOutputString (self)
            
            
             str = sprintf ('    drive : %d, %s,', self.label, self.type);
            
        end
        
    end
    
end