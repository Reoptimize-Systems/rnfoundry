classdef frictionModel < mbdyn.pre.base
    % base class for friction models
    
    properties
        frictionFunction;
        modelType;
    end
    
    methods
        
        function self = frictionModel (friction_function)
            
            
            assert (isa (friction_function, 'mbdyn.pre.scalarFunction'), ...
                'friction function must be an mbdyn.pre.scalarFunction object');
            
            self.frictionFunction = friction_function;
            
        end
        
    end
    
end