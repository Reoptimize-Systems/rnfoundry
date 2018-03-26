classdef revolutePin < mbdyn.pre.singleNodeJoint
    
    
    properties (GetAccess = public, SetAccess = protected)
        
        relativeOffset;
        relativeOffsetReference;
        nodeRelativeOrientation;
        nodeRelativeOrientationReference;
            
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
            %  'NodeRelativeOrientation' - orientation of the pin relative
            %    to the attached node. If both this and PinOrientation are
            %    supplied they must match.
            %
            %  'NodeOffsetReference' - string defining the reference frame
            %    in which the relative_offset (see above) is defined. By
            %    default this is 'global'. Other possibilities are 'node',
            %    to define a point in the reference frame of the node,
            %    and 'local', which is the same as 'node'.
            %
            %  'NodeRelativeOrientationReference' - string defining the reference frame
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
            options.NodeRelativeOrientationReference = 'global';
            
            options = parse_pv_pairs (options, varargin);
            
            % call the superclass constructor
            self = self@mbdyn.pre.singleNodeJoint (node);
            
            self.type = 'revolute pin';
            
            self.checkCartesianVector (relative_offset, true, 'relative_offset');
            self.checkCartesianVector (absolute_pin_position, true, 'absolute_pin_position');
            self.checkOrientationMatrix (options.NodeRelativeOrientation, true, 'NodeRelativeOrientation');
            self.checkOrientationMatrix (options.PinOrientation, true, 'PinOrientation');
            self.checkNodeReferenceType (options.NodeOffsetReference, true, 'NodeOffsetReference');
            self.checkNodeReferenceType (options.NodeRelativeOrientationReference, true, 'NodeRelativeOrientationReference');

            self.relativeOffset = relative_offset;
            self.relativeOffsetReference = options.NodeOffsetReference;
            
            self.nodeRelativeOrientation = options.NodeRelativeOrientation;
            self.nodeRelativeOrientationReference = options.NodeRelativeOrientationReference;

            self.pinPosition = absolute_pin_position;
            self.absolutePinOrientation = options.PinOrientation;
            
            self.initialTheta = options.InitialTheta;
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = generateMBDynInputString@mbdyn.pre.singleNodeJoint(self);
            
            str = self.addOutputLine (str, sprintf('%d', self.node.label), 2, true, 'node label');

            str = self.addOutputLine ( str, ...
                                       self.commaSepList ( 'position', ...
                                                           'reference', ...
                                                           self.relativeOffsetReference,...
                                                           self.relativeOffset ), ...
                                       3, ...
                                       true, ...
                                       'node relative position' );
            
            if ~isempty (self.nodeRelativeOrientation)
                
                str = self.addOutputLine ( str, ...
                                           self.commaSepList ( 'orientation', ...
                                                               'reference', ...
                                                               self.nodeRelativeOrientationReference, ...
                                                               self.nodeRelativeOrientation ), ...
                                           3, ...
                                           true, ...
                                           'node relative orientation' );
                                       
            end
            
            out = self.makeCellIfNot (self.pinPosition);
            addcomma = ~(isempty (self.absolutePinOrientation) && isempty (self.initialTheta));
            str = self.addOutputLine ( str, self.commaSepList ('position', out{:}), 2, addcomma, 'pin absolute position');
            
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
        
        function [ref_pos, ref_orient] = reference (self)
            % returns a reference object for the joint position
            
            switch self.relativeOffsetReference
                
                case {'node', 'local'}
                    posref = self.node.reference ();
                case {'global', ''}
                    posref = mbdyn.pre.globalref ();
                otherwise
                    error ('Unrecognised reference type');
                    
            end
            
            switch self.nodeRelativeOrientationReference
                
                case {'node', 'local'}
                    orientref = self.node.reference ();
                case  {'global', ''}
                    orientref = mbdyn.pre.globalref ();
                otherwise
                    error ('Unrecognised reference type');
                    
            end
            
            if ischar (self.relativeOffset)
                reloffset = [0;0;0];
            else
                reloffset = self.relativeOffset;
            end
            
            if ischar (self.nodeRelativeOrientation) ...
                || isempty (self.nodeRelativeOrientation)
                
                relorient = mbdyn.pre.orientmat ('eye');
                
            else
                relorient = mbdyn.pre.orientmat ('orientation', self.getOrientationMatrix (self.nodeRelativeOrientation));
            end
            
            ref_pos = mbdyn.pre.reference ( reloffset, ...
                                            relorient, ...
                                            [], ...
                                            [], ...
                                            'Parent', posref );
                                        
            ref_orient = mbdyn.pre.reference ( reloffset, ...
                                               relorient, ...
                                               [], ...
                                               [], ...
                                               'Parent', orientref );
                                    
        end
        
        function draw (self, varargin)
            
            options.AxesHandle = self.drawAxesH;
            options.ForceRedraw = false;
            options.Mode = 'solid';
            
            options = parse_pv_pairs (options, varargin);
            
%             draw@mbdyn.pre.element ( self, ...
%                 'AxesHandle', options.AxesHandle, ...
%                 'ForceRedraw', options.ForceRedraw, ...
%                 'Mode', options.Mode );
% 
%             self.setTransform ();

            if options.ForceRedraw
                self.needsRedraw = true;
            end
            
            self.checkAxes (options.AxesHandle);
            
            nodepos = self.node.absolutePosition;
            [jref_pos, ~] = self.reference ();
            jpos = jref_pos.pos ();
                
            if ~self.needsRedraw
                % always have to redraw line, can't just transform objects
                delete (self.shapeObjects{1})
                self.shapeObjects{1} =  line ( self.drawAxesH, ...
                                               [ nodepos(1), jpos(1) ], ...
                                               [ nodepos(2), jpos(2) ], ...
                                               [ nodepos(3), jpos(3) ], ...
                                               'Color', self.drawColour );
                                       
            end
            
            if isempty (self.shapeObjects) ...
                    || self.needsRedraw
                % a full redraw is needed (and not just a modification of
                % transform matrices for the objects).
                
                % delete the current patch object
                self.deleteAllDrawnObjects ();
                
                self.shapeObjects = { line( self.drawAxesH, ...
                                            [ nodepos(1), jpos(1) ], ...
                                            [ nodepos(2), jpos(2) ], ...
                                            [ nodepos(3), jpos(3) ], ...
                                            'Color', self.drawColour ), ...
                                      line( self.drawAxesH, ...
                                            [ 0, 0 ], ...
                                            [ 0, 0 ], ...
                                            [ -self.sz/2, self.sz/2 ], ...
                                            'Parent', self.transformObject, ...
                                            'Color', self.drawColour, ...
                                            'LineStyle', '--' )
                                    };
                
                self.needsRedraw = false;

%                 if options.Light
%                     light (self.drawAxesH);
%                 end
                
            end
            
            self.setTransform ();
            
        end
        
    end
    
    methods (Access = protected)
        
       function setTransform (self)
            
            [ref_pos, ref_orient] = self.reference ();
            
            M = [ ref_orient.orientm.orientationMatrix , ref_pos.pos; ...
                  0, 0, 0, 1 ];
            
            % matlab uses different convention to mbdyn for rotation
            % matrix
            M = self.mbdynOrient2Matlab (M);
                  
            set ( self.transformObject, 'Matrix', M );
            
        end 
        
    end
    
    
    
end