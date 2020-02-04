classdef drive < mbdyn.pre.base
    % generic base class for all mbdyn drives
    
    properties
        
    end
    
    methods
        
        function self = drive ()
            % construct a generic mbdyn drive object
            
            self = self@mbdyn.pre.base ();
            
        end
        
        function str = generateMBDynInputString (self)
            
             str = generateMBDynInputString@mbdyn.pre.element (self);
            
             str = sprintf ('%s    drive : %d, %s,', str, self.label, self.type);
            
        end
        
    end
    
end