classdef revoluteHinge < mbdyn.pre.twoNodeJoint
    
    
    properties (GetAccess = public, SetAccess = protected)
        
        relativeOffset1;
        relativeOrientation1;
        offset1Reference;
        orientation1Reference;
            
        relativeOffset2;
        relativeOrientation2;
        offset2Reference;
        orientation2Reference;
        
    end
    
    properties (Dependent)
        absoluteJointPosition;
    end
    
    methods
        
        function self = revoluteHinge (node1, node2, position1, position2, varargin)
            
            options.RelativeOrientation1 =  [];
            options.RelativeOrientation2 =  [];
            options.InitialTheta = [];
            options.Friction = [];
            options.Preload = [];
            options.FrictionModel = [];
            options.ShapeFunction = [];
            options.Offset1Reference = 'node';
            options.Offset2Reference = 'node';
            options.Orientation1Reference = 'node';
            options.Orientation2Reference = 'node';
            
            options = parse_pv_pairs (options, varargin);
            
            % call the superclass constructor
            self = self@mbdyn.pre.twoNodeJoint (node1, node2);
            
            self.type = 'revolute hinge';
            
            self.relativeOffset1 = self.checkJointPositionOffset ({options.Offset1Reference, position1});
            self.offset1Reference = options.Offset1Reference;
            
            self.relativeOrientation1 = self.checkJointOrientationOffset ({options.Orientation1Reference, options.RelativeOrientation1});
            self.orientation1Reference = options.Orientation1Reference;
            
            self.relativeOffset2 = self.checkJointPositionOffset ({options.Offset2Reference, position2});
            self.offset2Reference = options.Offset2Reference;
            
            self.relativeOrientation2 = self.checkJointOrientationOffset ({options.Orientation2Reference, options.RelativeOrientation2});
            self.orientation2Reference = options.Orientation2Reference;
            
        end
        
        function str = generateOutputString (self)
            
            str = generateOutputString@mbdyn.pre.twoNodeJoint(self);
            
            str = self.addOutputLine (str, sprintf('%d', self.node1.label), 2, true, 'node 1 label');
            
            out = self.makeCellIfNot (self.relativeOffset1);
            str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 3, true);
            
            if ~isempty (self.relativeOrientation1)
                out = self.makeCellIfNot (self.relativeOrientation1);
                str = self.addOutputLine (str, self.commaSepList ('orientation', out{:}), 3, true);
            end
            
            str = self.addOutputLine (str, sprintf('%d', self.node2.label), 2, true, 'node 2 label');
            
            out = self.makeCellIfNot (self.relativeOffset2);
            addcomma = ~isempty (self.relativeOrientation2);
            str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 3, addcomma);
            
            if ~isempty (self.relativeOrientation2)
                out = self.makeCellIfNot (self.relativeOrientation2);
                str = self.addOutputLine (str, self.commaSepList ('orientation', out{:}), 3, false);
            end
            
            str = self.addOutputLine (str, ';', 1, false, sprintf('end %s', self.type));
            
        end
        
        function draw (self, varargin)
            
            options.AxesHandle = [];
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
    
%     % gettters setters
%     methods
%         
%         function pos = get.absoluteJointPosition (self)
%             
%             if isa (self.offset1Reference, 'mbdyn.pre.reference')
%                 
%             elseif ischar (self.offset1Reference)
%                 
%                 
%             end
%             
%         end
%         
%     end
    
    methods (Access = protected)
        
        function setTransform (self)
            
            switch self.offset1Reference
                
                case ''
                    ref_pos_base = mbdyn.pre.reference (self.node1.absolutePosition, ...
                                                        self.node1.absoluteOrientation, ...
                                                        self.node1.absoluteVelocity, ...
                                                        self.node1.absoluteAngularVelocity);
                                                    
                case 'node'
                    ref_pos_base = mbdyn.pre.reference (self.node1.absolutePosition, ...
                                                        self.node1.absoluteOrientation, ...
                                                        self.node1.absoluteVelocity, ...
                                                        self.node1.absoluteAngularVelocity);
                                                    
                case 'local'
                    ref_pos_base = mbdyn.pre.reference (self.node1.absolutePosition, ...
                                                        self.node1.absoluteOrientation, ...
                                                        self.node1.absoluteVelocity, ...
                                                        self.node1.absoluteAngularVelocity);
                                                    
                case 'other node'
                    ref_pos_base = mbdyn.pre.reference (self.node2.absolutePosition, ...
                                                        self.node2.absoluteOrientation, ...
                                                        self.node2.absoluteVelocity, ...
                                                        self.node2.absoluteAngularVelocity);
                                                    
                case 'global'
                    ref_pos_base = mbdyn.pre.globalref;
            end
            
            switch self.orientation1Reference
                
                case ''
                    ref_orient_base = mbdyn.pre.reference (self.node1.absolutePosition, ...
                                                           self.node1.absoluteOrientation, ...
                                                           self.node1.absoluteVelocity, ...
                                                           self.node1.absoluteAngularVelocity);
                                                       
                case 'node'
                    ref_orient_base = mbdyn.pre.reference (self.node1.absolutePosition, ...
                                                           self.node1.absoluteOrientation, ...
                                                           self.node1.absoluteVelocity, ...
                                                           self.node1.absoluteAngularVelocity);
                case 'local'
                    ref_orient_base = mbdyn.pre.reference (self.node1.absolutePosition, ...
                                                           self.node1.absoluteOrientation, ...
                                                           self.node1.absoluteVelocity, ...
                                                           self.node1.absoluteAngularVelocity);
                                                       
                case 'other node'
                    ref_orient_base = mbdyn.pre.reference (self.node2.absolutePosition, ...
                                                           self.node2.absoluteOrientation, ...
                                                           self.node2.absoluteVelocity, ...
                                                           self.node2.absoluteAngularVelocity);
                                                       
                case 'global'
                    ref_orient_base = mbdyn.pre.globalref;
                    
            end

%                                         
            ref_joint_offset = mbdyn.pre.reference (self.relativeOffset1{end}, ...
                            mbdyn.pre.orientmat ('orientation', self.relativeOrientation1{end}), ...
                            [], ...
                            [], ...
                            'Parent', ref_pos_base);
            
            ref_joint_orient = mbdyn.pre.reference (self.relativeOffset1{end}, ...
                            mbdyn.pre.orientmat ('orientation', self.relativeOrientation1{end}), ...
                            [], ...
                            [], ...
                            'Parent', ref_orient_base);
                        
            M = [ ref_joint_orient.orientm.orientationMatrix, ref_joint_offset.pos; ...
                  0, 0, 0, 1 ];
            
            % matlab uses different convention to mbdyn for rotation
            % matrix
            M = self.mbdynOrient2Matlab (M);
                  
            set ( self.transformObject, 'Matrix', M );
            
        end
        
    end
    
end