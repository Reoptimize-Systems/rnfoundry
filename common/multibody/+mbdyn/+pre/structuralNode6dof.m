classdef structuralNode6dof < mbdyn.pre.structuralNode
    
    properties (GetAccess = public, SetAccess = protected)
        
        absoluteOrientation;
        orientationDescription;
        
        absoluteAngularVelocity;
        
        omegaRotates;
        
    end
    
    methods
        
        function self = structuralNode6dof (type, varargin)
            
            options.OrientationDescription = '';
            options.AbsolutePosition = [0;0;0];
            options.AbsoluteOrientation = eye (3);
            options.AbsoluteVelocity = [0;0;0];
            options.AbsoluteAngularVelocity = [0;0;0];
            options.Accelerations = [];
            
            options = parse_pv_pairs (options, varargin);
            
            switch type
                
                case 'static'
                    
                case 'dynamic'
                    
                case 'modal'
                    
                otherwise
                    error ('Unrecognised structural node type');
            end
            
            self = self@mbdyn.pre.structuralNode ( ...
                       'AbsoluteVelocity', options.AbsoluteVelocity, ...
                       'AbsolutePosition', options.AbsolutePosition, ...
                       'Accelerations', options.Accelerations);
                   
            self.type = type;
            
            self.absoluteOrientation = options.AbsoluteOrientation;
            
            self.absoluteAngularVelocity = options.AbsoluteAngularVelocity;
            
            if ~isempty (options.OrientationDescription)
                
                self.checkOrientationDescription (options.OrientationDescription, true);
            
            end
            
            self.orientationDescription = options.OrientationDescription;
            
        end
        
    end
    
    methods (Access = public)
        
        function str = generateOutputString (self)
            
            str = self.addOutputLine ('' , '', 1, false, '6 DOF structural node');
            
            % delete newline character and space from start
            str(1:2) = [];
            
            str = self.addOutputLine (str, sprintf('structural : %d, %s', self.label, self.type), 1, true, 'label, type');
            
            str = self.addOutputLine (str, self.commaSepList (self.absolutePosition), 2, true, 'absolute position');
            
            str = self.addOutputLine (str, self.commaSepList (self.absoluteOrientation), 2, true, 'absolute orientation');
            
            if ~isempty (self.orientationDescription)
                str = self.addOutputLine (str, self.commaSepList ('orientation description', self.orientationDescription), 3, true);
            end
            
            str = self.addOutputLine (str, self.commaSepList (self.absoluteVelocity), 2, true, 'absolute velocity');
            
            addcomma = ~isempty (self.accelerations);
            
            str = self.addOutputLine (str, self.commaSepList (self.absoluteAngularVelocity), 2, addcomma, 'absolute angular velocity');
            
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
    
end