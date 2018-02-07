classdef simplePlaneHingeShape < mbdyn.pre.shapeFunction
    
    properties
        radius;
    end
    
    methods
        
        function self = simplePlaneHingeShape (radius)
            
            self = self@mbdyn.pre.shapeFunction ('simple plane hinge');
            
            self.checkNumericScalar (radius, true, 'radius');
            
            self.radius = radius;

        end
        
        function str = generateOutputString (self)
            
            str = self.commaSepList ( self.fcnType, self.radius );
            
        end
        
    end
    
end