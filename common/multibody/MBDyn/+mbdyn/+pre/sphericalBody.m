classdef sphericalBody < mbdyn.pre.body
    
    properties
        
        
    end
    
    
    methods
        
        function self = sphericalBody (mass, R, node)
            
            mbdyn.pre.base.checkNumericScalar (mass, true, 'mass');
            mbdyn.pre.base.checkNumericScalar (R, true, 'R');
            mbdyn.pre.base.checkIsStructuralNode (node, true, 'node');
            
            cog = [0; 0; 0];
            
            Isphere = (2/5) * mass * R^2;
            
            self = self@mbdyn.pre.body (mass, cog, diag ([Isphere, Isphere, Isphere]), node);
            
            self.setSize (R, R, R);
    
        end
        
    end
    
end