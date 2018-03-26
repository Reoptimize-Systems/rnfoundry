classdef simpleShape < mbdyn.pre.shapeFunction
    % base class for friction models
    
    properties

    end
    
    methods
        
        function self = simpleShape ()
            
            self = self@mbdyn.pre.shapeFunction ('simple');
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = self.fcnType;              
            
        end
        
    end
    
end