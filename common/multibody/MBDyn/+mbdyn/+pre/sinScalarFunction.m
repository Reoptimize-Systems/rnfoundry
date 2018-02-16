classdef sinScalarFunction < mbdyn.pre.scalarFunction
    
    properties
        
        A;
        w;
        theta;
        
    end
    
    methods
        
        function self = sinScalarFunction (name, A, w, theta)
            % sin scalar function constructor
            %
            %
            % Syntax
            %
            % ssf = mbdyn.pre.sinScalarFunction (name, A, w)
            % ssf = mbdyn.pre.sinScalarFunction (name, A, w, theta)
            %
            % Description
            %
            % Implements a differentiable scalar function
            %
            %   y = A * sin (w * x + theta)
            %
            % Input
            %
            %  name - unique name for the function
            %
            %  A - amplitude of sin function
            %
            %  w - angular frequency of sin function
            %
            %  theta - (optional) phase of sin function, defaults to zero
            %    if not supplied.
            %
            % Output
            %
            %  ssf - mbdyn.pre.constScalarFunction object
            %
            %
            % See Also: mbdyn.pre.scalarFunction
            %
            
            if nargin < 4
                theta = 0;
            end
            
            self = self@mbdyn.pre.scalarFunction (name, 'sin');
            
            self.checkNumericScalar (A, true, 'A');
            self.checkNumericScalar (w, true, 'w');
            self.checkNumericScalar (theta, true, 'theta');
            
            self.A = A;
            self.w = w;
            self.theta = theta;
            
        end
        
        function str = generateOutputString (self)
            
            str = self.commaSepList ( ['"', self.name, '"'], ...
                                      self.fcnType, ...
                                      self.A, self.w, self.theta );
            
            
        end
        
    end
    
end