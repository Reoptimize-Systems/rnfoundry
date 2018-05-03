classdef subsystem < mbdyn.pre.base
    % subsystem: a collection of mbdyn elements representing a subsystem
    
    
    properties (Access=protected)
        
        structuralConnections;
        electricConnections;
        hydraulicConnections;
        
       references;
       data;
       problems;
       controlData;
       nodes;
       drivers;
       elements;
        
    end
    
    
    methods
        
        function self = subsystem ()
            
        end
        
        
        
        
    end
    
    
end