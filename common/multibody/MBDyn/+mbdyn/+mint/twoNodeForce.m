classdef twoNodeForce < mbdyn.mint.base
    
    properties
        
        referenceNode;
        otherNode;
        
    end
    
    methods
        
        function self = twoNodeForce (reference_node, other_node)
            % constructor for twoNodeForce base class
        
            if ~ ( isa (reference_node, 'mbdyn.pre.structuralNode') ...
                    || mbdyn.pre.base.checkScalarInteger (reference_node, false) )
                error ('reference_node must be an mbdyn.pre.structuralNode or an integer indicating the node label')
            end
            
            if ~ ( isa (other_node, 'mbdyn.pre.structuralNode') ...
                    || mbdyn.pre.base.checkScalarInteger (other_node, false) )
                error ('other_node must be an mbdyn.pre.structuralNode or an integer indicating the node label')
            end
            
            self.referenceNode = reference_node;
            self.otherNode = other_node;
        
        end
        
    end
    
end