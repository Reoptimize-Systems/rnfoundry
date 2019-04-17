classdef sphericalBody < mbdyn.pre.body
    
    properties
        
        
    end
    
    
    methods
        
        function self = sphericalBody (mass, R, node, varargin)
            
            options.Name = '';
            
            options = parse_pv_pairs (options, varargin);
            
            mbdyn.pre.base.checkNumericScalar (mass, true, 'mass');
            mbdyn.pre.base.checkNumericScalar (R, true, 'R');
            mbdyn.pre.base.checkIsStructuralNode (node, true, 'node');
            
            cog = [0; 0; 0];
            
            Isphere = (2/5) * mass * R^2;
            
            self = self@mbdyn.pre.body (mass, cog, diag ([Isphere, Isphere, Isphere]), node, 'Name', options.Name);
            
            self.setSize (R, R, R);
    
        end
        
    end
    
end