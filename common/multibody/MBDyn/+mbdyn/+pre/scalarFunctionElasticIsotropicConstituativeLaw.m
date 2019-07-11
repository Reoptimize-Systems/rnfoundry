classdef scalarFunctionElasticIsotropicConstituativeLaw < mbdyn.pre.constituativeLaw
    
    properties
        
        scalarFunction;
        
    end
    
    methods
        
        function self = scalarFunctionElasticIsotropicConstituativeLaw (scalar_function)
            % mbdyn.pre.scalarFunctionElasticIsotropicConstituativeLaw constructor
            %
            % Syntax
            %
            % sfeilaw = scalarFunctionElasticIsotropicConstituativeLaw (scalar_function)
            %
            % Description
            %
            % Applies stiffness according to a scalar function to represent
            % an analytical force-displacement curve of a single variable
            % that is automatically differentiated to compute the slope of
            % the curve, namely the local stiffness.
            %
            % Input
            %
            %  scalar_function - mbdyn.pre.scalarFunction object 
            %   representing the analytical force-displacement curve for
            %   this constituitive law.
            %
            % Output
            %
            %  sfeilaw - mbdyn.pre.scalarFunctionElasticIsotropicConstituativeLaw
            %   object
            %
            %
            % See Also:
            %
            
            assert ( isa (scalar_function, 'mbdyn.pre.scalarFunction'), ...
                     'scalar_function must be an mbdyn.pre.scalarFunction object' );
            
            self.type = 'scalar function elastic isotropic';
            self.scalarFunction = scalar_function;
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = self.commaSepList (self.type, self.scalarFunction.generateMBDynInputString ());
            
        end
        
    end
    
end