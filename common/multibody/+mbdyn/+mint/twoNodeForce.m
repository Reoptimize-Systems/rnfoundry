classdef twoNodeForce < mbdyn.mint.mintbase
    
    properties
        
        referenceNode;
        otherNode;
        
    end
    
    methods
        
        function self = twoNodeForce (reference_node, other_node)
        
            if ~ ( isa (reference_node, 'mbdyn.pre.structuralNode') ...
                    || isint2eps (reference_node) )
                error ('reference_node must be an mbdyn.pre.structuralNode or an integer indicating the node label')
            end
            
            if ~ ( isa (reference_node, 'mbdyn.pre.structuralNode') ...
                    || isint2eps (other_node) )
                error ('other_node must be an mbdyn.pre.structuralNode or an integer indicating the node label')
            end
            
            self.referenceNode = reference_node;
            self.otherNode = other_node;
        
        end
        
        function omat = referenceOrientation (self)
            
            
        end
        
    end
    
end