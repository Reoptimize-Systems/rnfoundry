classdef structuralNode3dof < mbdyn.pre.structuralNode
    
    properties
        
        
        
    end
    
    methods
        
        function self = structuralNode3dof (type, varargin)
            
            options.AbsolutePosition = [0;0;0];
            options.AbsoluteVelocity = [0;0;0];
            options.Accelerations = [];
            
            options = parse_pv_pairs (options, varargin);
            
            switch type
                
                case 'static displacement'
                    
                case 'dynamic displacement'
                    
                otherwise
                    error ('Unrecognised structural node type');
            end
            
            self = self@mbdyn.pre.structuralNode ( ...
                       'AbsoluteVelocity', options.AbsoluteVelocity, ...
                       'AbsolutePosition', options.AbsolutePosition, ...
                       'Accelerations', options.Accelerations);
                   
            self.type = type;
            
        end
        
        function str = generateOutputString (self)
            
            
            str = self.addOutputLine ('' , '', 1, false, '3 DOF structural node');
            
            % delete newline character and space from start
            str(1:2) = [];
            
            str = self.addOutputLine (str, sprintf('structural : %d, %s', self.label, self.type), 1, true, 'label, type');
            
            str = self.addOutputLine (str, self.commaSepList (self.absolutePosition), 2, true, 'absolute position');
            
            addcomma = ~isempty (self.accelerations);
            
            str = self.addOutputLine (str, self.commaSepList (self.absoluteVelocity), 2, addcomma, 'absolute velocity');
            
            if ~isempty (self.accelerations)
                
                if self.accelerations == true
                    str = self.addOutputLine (str, self.commaSepList ('accelerations', 'yes'), 2, false);
                else
                    str = self.addOutputLine (str, self.commaSepList ('accelerations', 'no'), 2, false);
                end
                
            end
            
            str = self.addOutputLine (str, ';', 1, false, 'end structural node');
            
        end
        
    end
    
    methods (Access = protected)
        function setTransform (self)
            
            M = [ eye(3), self.absolutePosition; ...
                  0, 0, 0, 1 ];
            
%             % matlab uses different convention to mbdyn for rotation
%             % matrix
%             M = self.mbdynOrient2Matlab (M);
%                   
            set ( self.transformObject, 'Matrix', M );
            
        end
    end
    
end