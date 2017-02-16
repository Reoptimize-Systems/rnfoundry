classdef singleNodeJoint < mbdyn.pre.joint
    
    
    properties (GetAccess = public, SetAccess = protected)
        
        node;
        
    end
    
    methods
        function self = singleNodeJoint (node)
        
            self.checkIsStructuralNode (node, true);
            
            self.node = node;
            
        end
    end
    
    methods (Access = protected)
        
        
    end
    
end