classdef followerCouple < mbdyn.pre.structuralCouple
    
    properties (GetAccess = public, SetAccess = protected)
        
        couple;
        
    end
    
    methods
        
        function self = followerCouple (node, couple_value, varargin)
            % <force_arglist> ::=
            % <node_label> ,
            % [ position , (Vec3) <relative_arm> , ]
            % (TplDriveCaller<Vec3>) <couple_value> 
            
            options.Position = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self = self@mbdyn.pre.structuralCouple (node, ...
                            'Position', options.Position);

            self.checkTplDriveCaller (couple_value, true, 'couple_value');
            
            self.couple = couple_value;
            
            self.type = 'follower';
            
        end
        
        function str = generateOutputString (self)
            
            str = generateOutputString@mbdyn.pre.couple (self);
            
            str = self.addOutputLine (str, sprintf('%d', self.node.label), 2, true);
            
            if ~isempty (self.position)
                str = self.addOutputLine (str, self.commaSepList ('position', self.position), 2, true);
            end
            
            str = self.addOutputLine (str, self.couple.generateOutputString(), 2, false);

            str = self.addOutputLine (str, ';', 1, false, 'end follower structural couple');
            
        end
        
    end
    
end