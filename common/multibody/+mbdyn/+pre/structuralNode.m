classdef structuralNode < mbdyn.pre.node
    
    properties (GetAccess = public, SetAccess = protected)
        
        absolutePosition;
        absoluteVelocity;
        
        % assembly
        positionInitialStiffness;
        velocityInitialStiffness;
        
        accelerations;
        
        type;
        
    end
    
    methods
        
        function self = structuralNode (varargin)
            
            options.AbsolutePosition = [0;0;0];
            options.AbsoluteVelocity = [0;0;0];
            options.Accelerations = [];
            
            options = parse_pv_pairs (options, varargin);
            
            if ~isempty (options.Accelerations)
                if islogical (options.Accelerations) && isscalar (options.Accelerations)
                    self.accelerations = options.Accelerations;
                else
                    error ('Accelerations should be a scalar boolean true/false');
                end
            end
            
            % TODO: add input checking code
            self.absolutePosition = options.AbsolutePosition;
            self.absoluteVelocity = options.AbsoluteVelocity;
            
            % TODO: find out what values of initial stiffness can be
            % supplied
            self.positionInitialStiffness = [];
            self.velocityInitialStiffness = [];
            
            
        end
        
    end
    
    methods (Access = protected)
        
    end
    
end

