classdef followerCouple < mbdyn.pre.structuralCouple
    
    properties (GetAccess = public, SetAccess = protected)
        
        couple;
        
    end
    
    methods
        
        function self = followerCouple (node, couple_value, varargin)
            % followerCouple constructor (a type of structural couple)
            %
            % Syntax
            %
            % fsc = mbdyn.pre.followerCouple (node, couple_value)
            % fsc = mbdyn.pre.followerCouple (..., 'Parameter', value)
            %
            % Description
            %
            % 
            %
            % Input
            %
            %  node - 
            %
            %  couple_value - 
            %
            % Output
            %
            %  fsc - 
            %
            %
            %
            % See Also: 
            %
            
            options.Position = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self = self@mbdyn.pre.structuralCouple (node, 'follower', couple_value, ...
                            'Position', options.Position);
            
        end
        
        function str = generateOutputString (self)
            
            str = generateOutputString@mbdyn.pre.structuralCouple (self);
            
        end
        
    end
    
end