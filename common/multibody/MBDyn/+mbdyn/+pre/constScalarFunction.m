classdef constScalarFunction < mbdyn.pre.scalarFunction
    
    properties (GetAccess = public, SetAccess = protected)
        c;
    end
    
    methods
        
        function self = constScalarFunction (name, c)
            % const scalar function constructor
            %
            %
            % Syntax
            %
            % csf = mbdyn.pre.constScalarFunction (name, c)
            %
            % Description
            %
            % Implements a scalar function with constant value.
            %
            % Input
            %
            %  name - unique name for the function
            %
            %  c - scalar value of the function
            %
            % Output
            %
            %  csf - mbdyn.pre.constScalarFunction object
            %
            %
            % See Also: mbdyn.pre.scalarFunction
            %
            
            
            assert (ischar (name), 'name should be a character vector');
            
            self = self@mbdyn.pre.scalarFunction (name, 'const');
            
            self.checkNumericScalar (c, true, 'c');
            
            self.name = name;
            self.c = c;
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = self.commaSepList ( ['"', self.name, '"'], ...
                                      self.fcnType, ...
                                      self.c );
                                     
            
        end
        
    end
    
end