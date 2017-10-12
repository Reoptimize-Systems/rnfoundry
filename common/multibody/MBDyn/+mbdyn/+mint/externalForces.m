classdef externalForces < handle
    % class for managing collections of external forces on mbdyn nodes
    % applied through the external structural force type
    %
    % Assists with the organsation of the nodes and forces 
    %
    %
    
    properties
        
        forces;
        nodes;
%         objMBCNodal;
        
    end
    
    methods
        
        function self = externalForces (varargin)
            
            options.PreProcExternalStructuralForce = [];
            options.PreProcSystem = [];
            
            options = parse_pv_pairs (options, varargin);
            
            
        end
        
        function forces = processForces (self)
            % calculates the forces on the external structural nodes at the
            % current time step
            
            
        end
        
    end
    
end