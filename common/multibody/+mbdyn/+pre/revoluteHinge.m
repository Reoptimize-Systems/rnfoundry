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
            
            self.checkJointPositionOffset (position1);
            self.checkJointPositionOffset (position2);
            self.checkJointOrientationOffset (options.RelativeOrientation1);
            self.checkJointOrientationOffset (options.RelativeOrientation2);
            
            self.checkNodeReferenceType (options.Offset1Reference, true);
            self.checkNodeReferenceType (options.Offset2Reference, true);
            self.checkNodeReferenceType (options.Orientation1Reference, true);
            self.checkNodeReferenceType (options.Orientation2Reference, true);
            
            self.offset1Reference = options.Offset1Reference;
            self.orientation1Reference = options.Orientation1Reference;
            self.offset2Reference = options.Offset2Reference;
            self.orientation2Reference = options.Orientation2Reference;
            
            self.relativeOffset1 = {'reference', self.offset1Reference, position1};
            self.relativeOrientation1 = {'reference', self.orientation1Reference, self.getOrientationMatrix(options.RelativeOrientation1)};

            self.relativeOffset2 = {'reference', self.offset2Reference, position2};
            self.relativeOrientation2 = {'reference', self.orientation2Reference, self.getOrientationMatrix(options.RelativeOrientation2)};
            
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
            
            str = self.addOutputLine (str, ';', 1, false, 'end revolute hinge');
            
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
                
                case 'node'
                    ref_pos_base = mbdyn.pre.reference (self.node1.absolutePosition, ...
                                                        self.node1.absoluteOrientation, ...
                                                        self.node1.absoluteVelocity, ...
                                                        self.node1.absoluteAngularVelocity);
                case 'global'
                    ref_pos_base = mbdyn.pre.globalref;
            end
            
            switch self.orientation1Reference
                
                case 'node'
                    ref_orient_base = mbdyn.pre.reference (self.node1.absolutePosition, ...
                                                           self.node1.absoluteOrientation, ...
                                                           self.node1.absoluteVelocity, ...
                                                           self.node1.absoluteAngularVelocity);
                case 'global'
                    ref_orient_base = mbdyn.pre.globalref;
                    
            end

%                                         
            ref_joint = mbdyn.pre.reference (self.relativeOffset1{end}, ...
                            mbdyn.pre.orientmat ('orientation', self.relativeOrientation1{end}), ...
                            [], ...
                            [], ...
                            'PositionParent', ref_pos_base, ...
                            'OrientParent', ref_orient_base);
            
            M = [ ref_joint.orientm.orientationMatrix , ref_joint.pos; ...
                  0, 0, 0, 1 ];
            
            % matlab uses different convention to mbdyn for rotation
            % matrix
            M = self.mbdynOrient2Matlab (M);
                  
            set ( self.transformObject, 'Matrix', M );
            
        end
        
    end
    
end