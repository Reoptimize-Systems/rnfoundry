classdef shapeFunction < mbdyn.pre.base
    % base class for friction models
    
    properties
        fcnType;
    end
    
    methods
        
        function self = shapeFunction (fcnType)
            
            assert (ischar (fcnType), ...
                'fcnType must be a chararray');
            
            self.fcnType = fcnType;
            
        end
        
    end
    
end