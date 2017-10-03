classdef linearDrive < mbdyn.pre.drive
    
    properties (GetAccess = public, SetAccess = private)
        constCoef;
        slopeCoef
    end
    
    methods
        
        function self = linearDrive (const_coef, slope_coef)
            
            self.checkNumericScalar (const_coef, true, 'const_coef');
            self.checkNumericScalar (slope_coef, true, 'slope_coef');
            
            self.constCoef = const_coef;
            self.slopeCoef = slope_coef;
            
            self.type = 'linear';
            
        end
        
        function str = generateOutputString (self)
            
            str = self.commaSepList (self.type, self.constCoef, self.slopeCoef);
            
        end
        
    end
    
end