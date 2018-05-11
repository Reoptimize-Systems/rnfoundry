classdef totalJoint < mbdyn.pre.twoNodeJoint
    
    properties (GetAccess = public, SetAccess = protected)
        
        relativeOffset1;
        relativeOffset1Reference;
        relativePositionOrientation1;
        relativePositionOrientation1Reference;
        relativeRotOrientation1;
        relativeRotOrientation1Reference;
            
        relativeOffset2;
        relativeOffset2Reference;
        relativePositionOrientation2;
        relativePositionOrientation2Reference;
        relativeRotOrientation2;
        relativeRotOrientation2Reference;
        
        positionStatus;
        orientationStatus;
        
    end
    
    
    methods
        
        function self = totalJoint (node1, node2, posstatus, orientstatus, varargin)
            
            options.RelativeOffset1 = [];
            options.RelativeOffset1Reference = 'node';
            options.RelativePositionOrientation1 =  [];
            options.RelativePositionOrientation1Reference = 'node';
            options.RelativeRotOrientation1 =  [];
            options.RelativeRotOrientation1Reference = 'node';
            
            options.RelativeOffset2 =  [];
            options.RelativeOffset2Reference = 'node';
            options.RelativePositionOrientation2 =  [];
            options.RelativePositionOrientation2Reference = 'node';
            options.RelativeRotOrientation2 =  [];
            options.RelativeRotOrientation2Reference = 'node';
            
            options = parse_pv_pairs (options, varargin);
            
            % call the superclass constructor
            self = self@mbdyn.pre.twoNodeJoint (node1, node2);
            
            self.type = 'total joint';
            
            if ~isempty (options.RelativeOffset1)
                self.relativeOffset1 = self.checkJointPositionOffset ({options.RelativeOffset1Reference, options.RelativeOffset1});
            else
                self.relativeOffset1 = [];
            end
            
            if ~isempty (options.RelativePositionOrientation1)
                self.relativePositionOrientation1 = self.checkJointOrientationOffset ({options.RelativePositionOrientation1Reference, options.RelativePositionOrientation1});
            else
                self.relativePositionOrientation1 = [];
            end
            
            if ~isempty (options.RelativeRotOrientation1)
                self.relativeRotOrientation1 = self.checkJointOrientationOffset ({options.RelativeRotOrientation1Reference, options.RelativeRotOrientation1});
            else
                self.relativeRotOrientation1 = [];
            end
            
            if ~isempty (options.RelativeOffset2)
                self.relativeOffset2 = self.checkJointPositionOffset ({options.RelativeOffset2Reference, options.RelativeOffset2});
            else
                self.relativeOffset2 = [];
            end
            
            if ~isempty (options.RelativePositionOrientation2)
                self.relativePositionOrientation2 = self.checkJointOrientationOffset ({options.RelativePositionOrientation2Reference, options.RelativePositionOrientation2});
            else
                self.relativePositionOrientation2 = [];
            end
            
            if ~isempty (options.RelativeRotOrientation2)
                self.relativeRotOrientation2 = self.checkJointOrientationOffset ({options.RelativeRotOrientation2Reference, options.RelativeRotOrientation2});
            else
                self.relativeRotOrientation2 = [];
            end
            
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
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for totalJoint
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (tj)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  tpj - mbdyn.pre.totalJoint object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = generateMBDynInputString@mbdyn.pre.twoNodeJoint (self);
            
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
            
            str = self.addOutputLine (str, ';', 1, false, sprintf ('end %s', self.type));
            
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