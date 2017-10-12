classdef structuralCouple < mbdyn.pre.couple
    
    properties (GetAccess = public, SetAccess = protected)
        
        position;
        postionReference;
        node;
        
    end
    
    methods
        
        function self = structuralCouple (node, varargin)
            
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
            
            str = generateOutputString@mbdyn.pre.couple (self);
            
        end
        
    end
    
end