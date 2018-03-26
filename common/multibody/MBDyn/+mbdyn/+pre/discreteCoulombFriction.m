classdef discreteCoulombFriction < mbdyn.pre.frictionModel
    
    properties
        
        sigma2;
        velocityRatio;
        
    end
    
    methods
        
        function self = discreteCoulombFriction (friction_fcn, varargin)
            % Coulomb model with viscous friction and internal states to resolve stick/slip conditions.
            %
            % Syntax
            %
            % self = discreteCoulombFriction (scalar_fcn, varargin)
            %
            % Description
            %
            % discreteCoulombFriction implements a Coulomb model with
            % viscous friction and internal states to resolve stick/slip
            % conditions.
            %
            % Note that discrete coulomb cannot handle stick conditions
            % when the reaction force in the constraint goes to zero
            % (resulting in a singular Jacobian matrix). The preload
            % parameter needs to be used to make sure the joint is
            % preloaded as appropriate (see the documentation of the
            % specific joint for details).
            %
            % Input
            %
            %  friction_fcn - scalar function object (any class derived from
            %    mbdyn.pre.scalarFunction)
            %
            % Additional options may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'Sigma2' - 
            %
            %  'VelocityRatio' - scalar value used to discriminate
            %    stick/slip conditions. Default is 0.8 if not supplied.
            %
            % Output
            %
            %
            %
            % See Also: 
            %

            options.Sigma2 = [];
            options.VelocityRatio = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self = self@mbdyn.pre.frictionModel (friction_fcn);
            
            if ~isempty (options.Sigma2) 
                self.checkNumericScalar (options.Sigma2, true, 'Sigma2');
            end
            if ~isempty (options.VelocityRatio) 
                self.checkNumericScalar (options.VelocityRatio, true, 'VelocityRatio');
            end
            
            self.modelType = 'discrete coulomb';
            self.sigma2 = options.Sigma2;
            self.velocityRatio = options.VelocityRatio;
            
        end
        
        function str = generateMBDynInputString (self)
            

            str = [ self.modelType, ',' ];
            
            addcomma = ~isempty (self.sigma2) ...
                || ~isempty (self.velocityRatio);
            
            str = self.addOutputLine ( str, ...
                                       self.frictionFunction.generateMBDynInputString (), ...
                                       1, addcomma ); % don't add a comment!
            
            addcomma = ~isempty (self.velocityRatio);
            
            if ~isempty (self.sigma2) 
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ('sigma2', self.sigma2), ...
                                           1, addcomma );
            end
            
            if ~isempty (self.velocityRatio) 
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ('velocity ratio', self.velocityRatio), ...
                                           1, false );
            end
            
        end
        
        
    end
    
end