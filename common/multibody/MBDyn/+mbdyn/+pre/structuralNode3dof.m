classdef structuralNode3dof < mbdyn.pre.structuralNode
    
    properties
        
        
        
    end
    
    methods
        
        function self = structuralNode3dof (type, varargin)
            
            options.AbsolutePosition = [0;0;0];
            options.AbsoluteVelocity = [0;0;0];
            options.Accelerations = [];
            options.HumanReadableLabel = '';
            options.Scale = [];
            options.Output = [];
            
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
                       'Accelerations', options.Accelerations, ...
                       'HumanReadableLabel', options.HumanReadableLabel, ...
                       'Scale', options.Scale, ...
                       'Output', options.Output);
                   
            self.type = type;
            
        end
        
        function str = generateOutputString (self)
            
            
            nodestr = generateOutputString@mbdyn.pre.structuralNode (self);
            
            str = self.addOutputLine ('' , '', 1, false, '3 DOF structural node');
            
            % delete newline character and space from start
            str(1:2) = [];
            
            str = self.addOutputLine (str, sprintf('structural : %d, %s', self.label, self.type), 1, true, 'label, type');
            
            str = self.addOutputLine (str, self.commaSepList ('position', self.absolutePosition), 2, true, 'absolute position');
            
            addcomma = ~isempty (self.accelerations) || ~isempty (nodestr);
            
            str = self.addOutputLine (str, self.commaSepList ('velocity', self.absoluteVelocity), 2, addcomma, 'absolute velocity');
            
            addcomma = ~isempty (nodestr);
            
            if ~isempty (self.accelerations)
                
                if self.accelerations == true
                    str = self.addOutputLine (str, self.commaSepList ('accelerations', 'yes'), 2, addcomma);
                else
                    str = self.addOutputLine (str, self.commaSepList ('accelerations', 'no'), 2, addcomma);
                end
                
            end
            
            if ~isempty (nodestr)
                str = self.addOutputLine (str, nodestr, 2, false);
            end
            
            str = self.addOutputLine (str, ';', 1, false, 'end structural node');
            
        end
        
        function ref = reference (self)
            % returns an mbdyn.pre.reference for the node
            
            ref = mbdyn.pre.reference ( self.absolutePosition, ...
                                        [], ...
                                        self.absoluteVelocity, ...
                                        [] );

        end
        
        function abspos = relativeToAbsolutePosition (self, pos)
            % convert a position in the reference frame of the node to global
            
            self.checkCartesianVector (pos);
            
            ref_node = reference (self);
                                         
            ref_out = mbdyn.pre.reference ( pos, ...
                                            [], ...
                                            [], ...
                                            [], ...
                                            ref_node );
                                        
            abspos = ref_out.position;
            
        end
        
        function absorientm = relativeToAbsoluteOrientation (self, orientation)
            % convert an orientation in the reference frame of the node to global
            
            self.checkOrientationMatrix (orientation);
            
            ref_node = reference (self);
                                         
            ref_out = mbdyn.pre.reference ( [], ...
                                            orientation, ...
                                            [], ...
                                            [], ...
                                            ref_node );
                                        
            absorientm = ref_out.orientm;
            
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