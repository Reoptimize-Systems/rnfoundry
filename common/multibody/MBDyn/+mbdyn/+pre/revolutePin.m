classdef revolutePin < mbdyn.pre.singleNodeJoint
    
    
    properties (GetAccess = public, SetAccess = protected)
        
        relativeOffset;
        nodeRelativeOrientation;
            
        pinPosition;
        absolutePinOrientation;
        
        initialTheta;
    end
    
    
    methods
        
        function self = revolutePin (node, relative_offset, absolute_pin_position, varargin)
            % joint which only allows the absolute rotation of a node about
            % a given axis, which is axis 3 in the reference systems
            % defined by the two orientation statements
            %
            % Syntax
            %
            % obj = mbdyn.pre.revolutePin (node, relative_offset, absolute_pin_position)
            % obj = mbdyn.pre.revolutePin (..., 'Parameter', value)
            %
            % Description
            %
            % revolutePin is a joint which only allows the absolute
            % rotation of a node about a given axis, which is axis 3 in the
            % reference systems defined by the two orientation statements,
            % i.e. rotation occurs about the local axis 3 of the
            % constraint.
            %
            % Input
            %
            %  node - mbdyn.pre.structuralNode object (or child class).
            %    This is the node which will be constrained by the pin
            %    joint.
            %
            %  relative_offset - the position of the node relative to the
            %    hinge point. Note that by default this is a position
            %    relative to the global coordinate system, i.e. it is
            %    actually just the absolute position of the hinge point.
            %    You can specify positions relative to other coordinate
            %    systems using the 'NodeOffsetReference' option described
            %    below. 
            %
            %  absolute_pin_position - the position of the hinge point in
            %    the global reference frame.
            %
            % Additional optional arguments can be supplied using
            % parameter-value pairs. The available options are:
            %
            %  'PinOrientation' - orientation of the pin joint. Rotation
            %    occurs about the local axis 3 of this orientation. Can be
            %    supplied as an mbdyn.pre.orientmat object or a 3x3
            %    orientation matrix.
            %
            %  'NodeRelativeOrientation' - 
            %
            %  'NodeOffsetReference' - string defining the reference frame
            %    in which the relative_offset (see above) is defined. By
            %    default this is 'global'. Other possibilities are 'node',
            %    to define a point in the reference frame of the node,
            %    and 'local', which is the same as 'node'.
            %
            %  'NodeOrientationReference' - string defining the reference frame
            %    in which the orientation of the joint relative to the node
            %    (see NodeRelativeOrientation above) is defined. By default
            %    this is 'global'. Other possibilities are 'node', to
            %    define a point in the reference frame of the node, and
            %    'local', which is the same as 'node'.
            %
            %  'InitialTheta' - The initial angular displacement of the
            %    joint.
            %
            % Output
            %
            %  obj - mbdyn.pre.revolutePin object
            %
            % See Also: 
            %
            % mbdyn.pre.revoluteHinge
            %
            %
            
            options.PinOrientation =  [];
            options.NodeRelativeOrientation = [];
            options.InitialTheta = [];
            options.NodeOffsetReference = 'global';
            options.NodeOrientationReference = 'global';
            
            options = parse_pv_pairs (options, varargin);
            
            % call the superclass constructor
            self = self@mbdyn.pre.singleNodeJoint (node);
            
            self.type = 'revolute pin';
            
            self.checkCartesianVector (relative_offset, true);
            self.checkCartesianVector (absolute_pin_position, true);
            self.checkOrientationMatrix (options.NodeRelativeOrientation, true);
            self.checkOrientationMatrix (options.PinOrientation, true);
            self.checkNodeReferenceType (options.NodeOffsetReference, true);
            self.checkNodeReferenceType (options.NodeOrientationReference, true);

            self.relativeOffset = {'reference', options.NodeOffsetReference, relative_offset};
            self.nodeRelativeOrientation = {'reference', options.NodeOrientationReference, self.getOrientationMatrix(options.NodeRelativeOrientation)};

            self.pinPosition = absolute_pin_position;
            self.absolutePinOrientation = self.getOrientationMatrix (options.PinOrientation);
            
            self.initialTheta = options.InitialTheta;
            
        end
        
        function str = generateOutputString (self)
            
            str = generateOutputString@mbdyn.pre.singleNodeJoint(self);
            
            str = self.addOutputLine (str, sprintf('%d', self.node.label), 2, true, 'node label');
            
            out = self.makeCellIfNot (self.relativeOffset);
            str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 3, true, 'node relative position' );
            
            if ~isempty (self.nodeRelativeOrientation)
                out = self.makeCellIfNot (self.nodeRelativeOrientation);
                str = self.addOutputLine (str, self.commaSepList ('orientation', out{:}), 3, true, 'node relative orientation');
            end
            
            out = self.makeCellIfNot (self.pinPosition);
            addcomma = ~(isempty (self.absolutePinOrientation) && isempty (self.initialTheta));
            str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 2, addcomma, 'pin absolute position');
            
            if ~isempty (self.absolutePinOrientation)
                addcomma = ~isempty (self.initialTheta);
                out = self.makeCellIfNot (self.absolutePinOrientation);
                str = self.addOutputLine (str, self.commaSepList ('orientation', out{:}), 2, addcomma, 'pin absolute orientation');
            end
            
            if ~isempty (self.initialTheta)
                out = self.makeCellIfNot (self.initialTheta);
                str = self.addOutputLine (str, self.commaSepList ('initial theta', out{:}), 2, false, 'initial theta');
            end
            
            str = self.addOutputLine (str, ';', 1, false, sprintf('end %s', self.type));
            
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
            
            ref_pin = mbdyn.pre.reference ( self.pinPosition, ...
                                            mbdyn.pre.orientmat ('orientation', self.absolutePinOrientation), ...
                                            [], []);
                                        
%             ref_joint = mbdyn.pre.reference (self.relativeOffset, ...
%                 mbdyn.pre.orientmat ('orientation', self.nodeRelativeOrientation), ...
%                 [], [], 'Parent', ref_pin);
            
            M = [ ref_pin.orientm.orientationMatrix , ref_pin.pos; ...
                  0, 0, 0, 1 ];
            
            % matlab uses different convention to mbdyn for rotation
            % matrix
            M = self.mbdynOrient2Matlab (M);
                  
            set ( self.transformObject, 'Matrix', M );
            
        end
        
    end
    
    
    
end