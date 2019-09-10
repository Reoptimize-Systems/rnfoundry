classdef sumScalarFunction < mbdyn.pre.scalarFunction
    
    properties (GetAccess = public, SetAccess = protected)
        f1;
        f2;
    end
    
    methods
        
        function self = sumScalarFunction (name, f1, f2)
            % const scalar function constructor
            %
            %
            % Syntax
            %
            % ssf = mbdyn.pre.sumScalarFunction (name, c)
            %
            % Description
            %
            % Implements a scalar function which is the sum of two other
            % scalar functions.
            %
            % Input
            %
            %  name - unique name for the function
            %
            %  f1 - mbdyn.pre.scalarFunction derived object, first scalar
            %   function to be summed
            %
            %  f2 - mbdyn.pre.scalarFunction derived object, second scalar 
            %   function to be summed
            %
            % Output
            %
            %  ssf - mbdyn.pre.sumScalarFunction object
            %
            %
            % See Also: mbdyn.pre.scalarFunction
            %
            
            
            assert (ischar (name), 'name should be a char array');
            
            self = self@mbdyn.pre.scalarFunction (name, 'const');
            
            assert ( isa (f1, 'mbdyn.pre.scalarFunction'), ...
                     'f1 must be a mbdyn.pre.scalarFunction derived object');
            assert ( isa (f2, 'mbdyn.pre.scalarFunction'), ...
                     'f2 must be a mbdyn.pre.scalarFunction derived object');
                 
            self.name = name;
            self.fcnType = 'sum';
            self.f1 = f1;
            self.f2 = f2;
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = self.commaSepList ( ['"', self.name, '"'], ...
                                      self.fcnType, ...
                                      self.f1.generateMBDynInputString (),  ...
                                      self.f2.generateMBDynInputString () );
                                     
            
        end
        
    end
    
end