classdef structuralForce < mbdyn.pre.force
    
    properties (GetAccess = public, SetAccess = protected)
        
        position;
        postionReference;
        node;
        
    end
    
    methods
        
        function self = structuralForce (node, varargin)
            
            % 'Postion' - vector defining the offset with respect to the
            %   node of the point where the force is applied.
            
            options.Position = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkIsStructuralNode (node, true);
            
            self.checkCartesianVector (options.Position);
            
            self.node = node;
            self.position = options.Position;
            
        end
        
        function str = generateOutputString (self)
            
            str = generateOutputString@mbdyn.pre.force(self);
            
        end
        
    end
    
end