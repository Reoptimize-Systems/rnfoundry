classdef powertakeoff < handle
   
    properties
        
    end
    
    methods (Abstract)
        
        forceAndTorque (self);
        forceSize (self);
        initDataStructure (self);
        loggingInfo(self);
        
    end
    
end