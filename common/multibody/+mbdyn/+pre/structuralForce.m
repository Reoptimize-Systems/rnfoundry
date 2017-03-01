classdef structuralForce < mbdyn.pre.force
    
    properties (GetAccess = public, SetAccess = protected)
        
    end
    
    methods
        
        function str = generateOutputString (self)
            
            str = generateOutputString@mbdyn.pre.force(self);
            
        end
        
    end
    
end