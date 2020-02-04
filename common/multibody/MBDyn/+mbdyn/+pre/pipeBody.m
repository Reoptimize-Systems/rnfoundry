classdef pipeBody < mbdyn.pre.body
    
    properties
        
        
    end
    
    
    methods
        
        function self = pipeBody (mass, R, h, aligned_axis, node, varargin)
            
            options.Name = '';
            
            options = parse_pv_pairs (options, varargin);
            
            mbdyn.pre.base.checkNumericScalar (mass, true, 'mass');
            mbdyn.pre.base.checkNumericScalar (R, true, 'R');
            mbdyn.pre.base.checkNumericScalar (h, true, 'h');
            mbdyn.pre.base.checkScalarInteger (aligned_axis, true, 'aligned_axis');
            mbdyn.pre.base.checkIsStructuralNode (node, true, 'node');
            
            cog = [0; 0; 0];

            Iax = 0.5 * mass .* R^2;
            
            Iother = (1/12) * mass * h^2 + 0.25 * mass .* R^2;
            
            switch aligned_axis
                
                case 1
                    
                    inertia_mat = diag ([Iax, Iother, Iother]);
                    
                case 2
                    
                    inertia_mat = diag ([Iother, Iax, Iother]);
                    
                case 3
                    
                    inertia_mat = diag ([Iother, Iother, Iax]);
                    
                otherwise
                    
                    error ('aligned_axis must be 1, 2 or 3');
                        
            end
            
            self = self@mbdyn.pre.body (mass, cog, inertia_mat, node, 'DefaultShape', 'cylinder', 'Name', options.Name);
            
            self.setSize (R, h);
    
        end
        
    end
    
end