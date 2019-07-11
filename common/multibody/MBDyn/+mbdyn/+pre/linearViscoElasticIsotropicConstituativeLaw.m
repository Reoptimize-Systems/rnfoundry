classdef linearViscoElasticIsotropicConstituativeLaw < mbdyn.pre.constituativeLaw
    
    properties
        
        stiffness;
        viscosity;
        propFactor;
        
    end
    
    methods
        
        function self = linearViscoElasticIsotropicConstituativeLaw (stiffness, viscosity, prop, varargin)
            % constructor for linearViscoElasticIsotropicConstituativeLaw
            %
            % Syntax
            %
            % lveilaw = linearViscoElasticIsotropicConstituativeLaw (stiffness, viscosity)
            % lveilaw = linearViscoElasticIsotropicConstituativeLaw (stiffness, viscosity, prop)
            %
            % Description
            %
            % Applies linear, istropic stiffness and damping coefficients.
            % The viscosity can be specified as a fixed value, or as a
            % factor of the stiffness.
            %
            % Input
            %
            %  stiffness - scalar value of the stiffness coefficient
            %
            %  viscosity - scalar value of the viscosity coefficient, or
            %   the keyword 'proportional'. If it is 'proportional', the
            %   value of the proportionality factor must be supplied in
            %   prop.
            %
            % Output
            %
            %  lveilaw -
            %
            %
            % See Also:
            %
            
            [options, nopass_list] = mbdyn.pre.linearViscoElasticIsotropicConstituativeLaw.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            self = self@mbdyn.pre.constituativeLaw (pvpairs{:});
            
            
            self.checkNumericScalar (stiffness, true, 'stiffness');
            
            if ischar (viscosity)
                if strcmp (viscosity, 'proportional')
                    if nargin > 2
                        self.checkNumericScalar (prop, false, 'prop');
                        self.propFactor = prop;
                        self.viscosity = [];
                    else
                        error ('viscosity was the keyword ''proportional'', but no value was suppled in prop');
                    end
                else
                    error ('viscosity must be a scalar value or the keyword ''proportional''');
                end
            elseif self.checkNumericScalar (viscosity, false, 'viscosity')
                self.viscosity = viscosity;
                self.propFactor = [];
            else
                error ('viscosity must be a scalar value or the keyword ''proportional''');
            end
            
            self.type = 'linear viscoelastic isotropic';
            self.stiffness = stiffness;
            
        end
        
        function str = generateMBDynInputString (self)
            
            if isempty (self.propFactor)
                specific_law_str = self.commaSepList (self.type, self.stiffness, self.viscosity);
            else
                specific_law_str = self.commaSepList (self.type, self.stiffness, 'proportional', self.propFactor);
            end
            
            str = generateMBDynInputString@mbdyn.pre.constituativeLaw (self, specific_law_str);
            
        end
        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.constituativeLaw.defaultConstructorOptions ();
            
            nopass_list = {};
            
        end
        
    end
    
end