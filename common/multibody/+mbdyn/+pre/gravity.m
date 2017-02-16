classdef gravity < mbdyn.pre.element
    
    properties
        gravityAcceleration;
        gravityModel;
%         origin;
    end
    
    methods
        
        function self = gravity (varargin)
            
            options.GravityModel = 'uniform';
            options.GravityAcceleration = mbdyn.pre.const (9.81);
            options.Origin = [];
            
            options = parse_pv_pairs (options, varargin);
            
            
            self.gravityAcceleration = options.GravityAcceleration;

        end
        
    end
    
end