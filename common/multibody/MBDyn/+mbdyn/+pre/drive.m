classdef drive < mbdyn.pre.base
    % generic base class for all mbdyn drives
    
    properties
        
    end
    
    methods
        
        function self = drive ()
            % construct a generic mbdyn drive object
            
            
        end
        
        function str = generateMBDynInputString (self)
            
             str = sprintf ('    drive : %d, %s,', self.label, self.type);
            
        end
        
    end
    
end