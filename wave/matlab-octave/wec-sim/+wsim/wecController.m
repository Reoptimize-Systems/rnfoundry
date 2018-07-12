classdef wecController < handle
    
    properties
        
        wecSimObj; % the wsim.wecSim object running the simulation
        logger; % wsim.logger object
        
    end
    
    
    methods
        
        function value = ptoControlOutput (self, id, pto_internal_variables)
            
            value = 0;
            
        end
        
        function start (self)
            
            
        end
        
        function advanceStep (self)
            
            
        end
        
        function finish (self)
            
            
        end
        
    end
    
end