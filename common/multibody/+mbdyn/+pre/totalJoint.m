classdef totalJoint < mbdyn.pre.twoNodeJoint
    
    properties (GetAccess = public, SetAccess = protected)
        
        relativeOffset1;
        relativePositionOrientation1;
        relativeRotOrientation1;
            
        relativeOffset2;
        relativePositionOrientation2;
        relativeRotOrientation2;
        
        positionStatus;
        orientationStatus;
        
    end
    
    
    methods
        
        function self = totalJoint (node1, node2, posstatus, orientstatus, varargin)
            
            options.RelativeOffset1 = [];
            options.RelativePositionOrientation1 =  [];
            options.RelativeRotOrientation1 =  [];
            
            options.RelativeOffset2 =  [];
            options.RelativePositionOrientation2 =  [];
            options.RelativeRotOrientation2 =  [];
            
            options = parse_pv_pairs (options, varargin);
            
            % call the superclass constructor
            self = self@mbdyn.pre.twoNodeJoint (node1, node2);
            
            self.type = 'total joint';
            
            self.checkJointPositionOffset (options.RelativeOffset1);
            self.checkJointPositionOffset (options.RelativeOffset2);
            self.checkJointOrientationOffset (options.RelativePositionOrientation1);
            self.checkJointOrientationOffset (options.RelativePositionOrientation2);
            
            self.relativeOffset1 = options.RelativeOffset1;
            self.relativePositionOrientation1 = self.getOrientationMatrix (options.RelativePositionOrientation1);
            self.relativeRotOrientation1 = self.getOrientationMatrix (options.RelativeRotOrientation1);

            self.relativeOffset2 = options.RelativeOffset2;
            self.relativePositionOrientation2 = self.getOrientationMatrix (options.RelativePositionOrientation2);
            self.relativeRotOrientation2 = self.getOrientationMatrix (options.RelativeRotOrientation2);
            
            self.checkPosStatus (posstatus);
            self.checkOrientStatus (orientstatus);
            
            if islogical (posstatus)
                if posstatus == true
                    posstatus = 'active';
                else
                    posstatus = 'inactive';
                end
            end
            
            self.positionStatus = posstatus;
            
            if islogical (orientstatus)
                if orientstatus == true
                    orientstatus = 'active';
                else
                    orientstatus = 'inactive';
                end
            end
            
            self.orientationStatus = orientstatus;
            
        end
        
        function str = generateOutputString (self)
            
            str = generateOutputString@mbdyn.pre.twoNodeJoint (self);
            
            str = self.addOutputLine (str, sprintf('%d', self.node1.label), 2, true, 'node 1 label');
            
            if ~isempty (self.relativeOffset1)
                out = self.makeCellIfNot (self.relativeOffset1);
                str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 3, true);
            end
            
            if ~isempty (self.relativePositionOrientation1)
                out = self.makeCellIfNot (self.relativePositionOrientation1);
                str = self.addOutputLine (str, self.commaSepList ('position orientation', out{:}), 3, true);
            end
            
            if ~isempty (self.relativeRotOrientation1)
                out = self.makeCellIfNot (self.relativeRotOrientation1);
                str = self.addOutputLine (str, self.commaSepList ('rotation orientation', out{:}), 3, true);
            end
            
            str = self.addOutputLine (str, sprintf('%d', self.node2.label), 2, true, 'node 2 label');
            
            if ~isempty (self.relativeOffset2)
                out = self.makeCellIfNot (self.relativeOffset2);
                str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 3, true);
            end
            
            if ~isempty (self.relativePositionOrientation2)
                out = self.makeCellIfNot (self.relativePositionOrientation2);
                str = self.addOutputLine (str, self.commaSepList ('position orientation', out{:}), 3, true);
            end
            
            if ~isempty (self.relativeRotOrientation2)
                out = self.makeCellIfNot (self.relativeRotOrientation2);
                str = self.addOutputLine (str, self.commaSepList ('rotation orientation', out{:}), 3, true);
            end
            
            str = self.addOutputLine (str, self.positionStatus, 2, true, 'position constraint status');
            
            str = self.addOutputLine (str, self.orientationStatus, 2, false, 'orientation constraint status');
            
            str = self.addOutputLine (str, ';', 1, false, 'end total joint');
            
        end
        
        function draw (self, varargin)
            
            options.AxesHandle = self.drawAxesH;
            options.ForceRedraw = false;
            options.Mode = 'solid';
            
            options = parse_pv_pairs (options, varargin);
            
            draw@mbdyn.pre.element ( self, ...
                'AxesHandle', options.AxesHandle, ...
                'ForceRedraw', options.ForceRedraw, ...
                'Mode', options.Mode );

            self.setTransform ();
            
        end
        
    end
    
    methods (Access = protected)
        
        function setTransform (self)
            
            ref_node = mbdyn.pre.reference (self.node1.absolutePosition, ...
                                            self.node1.absoluteOrientation, ...
                                            [], []);
                                        
            ref_joint = mbdyn.pre.reference (self.relativeOffset1, self.relativePositionOrientation1, [], [], 'Parent', ref_node);
            
            M = [ ref_joint.orientm.orientationMatrix , ref_joint.pos; ...
                  0, 0, 0, 1 ];
            
            % matlab uses different convention to mbdyn for rotation
            % matrix
            M = self.mbdynOrient2Matlab (M);
                  
            set ( self.transformObject, 'Matrix', M );
            
        end
        
        function checkPosStatus (self, posstatus)
            % checks if the position status choice is valid
            
            if (ischar (posstatus) && ~any(strcmp (posstatus, {'inactive', 'active', 'position', 'velocity'}))) ...
                    && ~islogical (posstatus)
                
                error ('total joint position status must be: inactive | active | position | velocity | boolean');
            end
            
        end
        
        function checkOrientStatus (self, orientstatus)
            % checks if the orientation status choice is valid
            
            if (ischar (orientstatus) && ~any(strcmp (orientstatus, {'inactive', 'active', 'rotation', 'angular velocity'}))) ...
                    && ~islogical (orientstatus)
                
                error ('total joint orientation status must be: inactive | active | rotation | angular velocity | boolean');
            end
            
        end
        
    end
    
    
    
end