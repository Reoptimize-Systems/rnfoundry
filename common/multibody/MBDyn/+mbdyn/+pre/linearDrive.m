classdef linearDrive < mbdyn.pre.drive
    
    properties (GetAccess = public, SetAccess = private)
        constCoef;
        slopeCoef
    end
    
    methods
        
        function self = linearDrive (const_coef, slope_coef)
            % evaluates the function y = mx + c on caller value
            %
            % Syntax
            %
            % ld = linearDrive (const_coef, slope_coef)
            %
            % Description
            %
            % Linear drive applies the function y = mx + c to the value
            % supplied by it's caller.
            %
            % Input
            %
            %  const_coef - constant coefficient of function
            %
            %  slope_coef - slope coefficient of function
            %
            % Output
            %
            %  ld - mbdyn.pre.linearDrive object
            %
            % See Also: 
            %

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