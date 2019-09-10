classdef scalarFunctionElasticIsotropicConstituativeLaw < mbdyn.pre.constituativeLaw
    
    properties
        
        scalarFunction;
        
    end
    
    methods
        
        function self = scalarFunctionElasticIsotropicConstituativeLaw (scalar_function, varargin)
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
            
            [options, nopass_list] = mbdyn.pre.scalarFunctionElasticIsotropicConstituativeLaw.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            assert ( isa (scalar_function, 'mbdyn.pre.scalarFunction'), ...
                     'scalar_function must be an mbdyn.pre.scalarFunction object' );

            % call superclass constructor
            self = self@mbdyn.pre.constituativeLaw ( pvpairs{:} );
            
            self.type = 'scalar function elastic isotropic';
            self.scalarFunction = scalar_function;
            
        end
        
        function str = generateMBDynInputString (self)

            str = generateMBDynInputString@mbdyn.pre.constituativeLaw (self, self.commaSepList (self.type, self.scalarFunction.generateMBDynInputString ()));
            
        end
        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.constituativeLaw.defaultConstructorOptions ();
            
            parentfnames = fieldnames (options);
            
            % would add new options here
            % options.NewOption = [];
            
            allfnames = fieldnames (options);
            
            nopass_list = setdiff (allfnames, parentfnames);
            
        end
        
    end
    
end