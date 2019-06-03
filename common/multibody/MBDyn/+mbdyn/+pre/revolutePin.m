classdef revolutePin < mbdyn.pre.singleNodeJoint
% class representing a revolute pin constraint (hinge pinned to global frame)
%
% Syntax
%
% obj = mbdyn.pre.revolutePin (node, relative_offset, absolute_pin_position)
% obj = mbdyn.pre.revolutePin (..., 'Parameter', value)
%
% Description
%
% revolutePin is a joint which only allows the absolute rotation of a node
% about a given axis, which is axis 3 in the reference systems defined by
% the two orientation statements, i.e. rotation occurs about the local axis
% 3 of the constraint.
%
% mbdyn.pre.revolutePin Methods:
%
%   revolutePin - mbdyn.pre.revolutePin constructor
%   defaultConstructorOptions - get the parent class's default options
%   draw - draw an mbdyn.pre.revolutePin in a figure
%   generateMBDynInputString - generates MBDyn input string for revolutePin object
%   reference - returns a mbdyn.pre.reference object for the joint position
%   setSize - set the size of the revolutePin in plots
%
    
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
            % mbdyn.pre.revolutePin constructor
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
            %  relative_offset - (3 x 1) vector (or 'null') representing the 
            %    position of the hinge point relative to the node. Note
            %    that by default this is a position relative to the global
            %    coordinate system, i.e. it is actually just the absolute
            %    position of the hinge. You can specify positions relative
            %    to other coordinate systems using the
            %    'NodeOffsetReference' option described below.
            %
            %  absolute_pin_position - (3 x 1) vector (or 'null') 
            %    representing the position of the hinge point in the global
            %    reference frame. This must match the position defined in
            %    'relative_offset'.
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
            
            [ options, nopass_list ] = mbdyn.pre.revolutePin.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs ( options, nopass_list);
            
            % call the superclass constructor
            self = self@mbdyn.pre.singleNodeJoint (node, 'DefaultShape', 'none', pvpairs{:});
            
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
            
            self.setSize (1);
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for revolutePin object
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (rp)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  rp - mbdyn.pre.revolutePin object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
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
            % returns a mbdyn.pre.reference object for the joint position
            
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
        
        function setSize (self, varargin)
            % set the size of the revolutePin in plots
            %
            % Syntax
            %
            % setSize (obj, radius)
            %
            % Description
            %
            % setSize is used to set the size of the default revolutePin
            % shape for plotting the revolutePin in a figure. This is used
            % when no STL file is available. The revolute pi in represented
            % as a dashed line parallel to the pin rotation axis passing
            % through the pin joint location, and a line linking the pin
            % joint location and attached node.
            %
            % Input
            %
            %  obj - mbdyn.pre.revolutePin object
            %
            %  pinlinelen - length of dashed pin axis line.
            %
            %
            % See Also: 
            %

            if self.stlLoaded
                
                setSize@mbdyn.pre.singleNodeJoint (varargin{:});
                
            else
                % cuboid, 3 arguments expected, x, y and z dimensions
                assert ( numel (varargin) == 1, ...
                         [ 'setSize requires one size input arguments when ', ...
                           'the shape is not from an STL file, the pin line length, ', ...
                           'which represent the bounding box of the shape' ] );

                self.checkNumericScalar (varargin{1}, true, 'pin line length');

                assert (varargin{1} > 0, 'pin line length must be greater than zero');

                self.shapeParameters(1) = varargin{1};
            end
        
            
        end
        
        
        function draw (self, varargin)
            % draw an mbdyn.pre.revolutePin in a figure
            
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
                                            [ -self.shapeParameters(1)/2, self.shapeParameters(1)/2 ], ...
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
            
%             % matlab uses different convention to mbdyn for rotation
%             % matrix
%             M = self.mbdynOrient2Matlab (M);
%                   
            set ( self.transformObject, 'Matrix', M );
            
        end 
        
    end
    
    methods (Static)
        
        function [ options, nopass_list ] = defaultConstructorOptions ()
            
            % get the parent class's default options 
            options = mbdyn.pre.singleNodeJoint.defaultConstructorOptions ();
            
            parentfnames = fieldnames (options);
            
            % add default options for revolutePin
            options.PinOrientation =  [];
            options.NodeRelativeOrientation = [];
            options.InitialTheta = [];
            options.NodeOffsetReference = 'global';
            options.NodeRelativeOrientationReference = 'global';
            
            allfnames = fieldnames (options);
            
            % get just the new option names
            C = setdiff (allfnames, parentfnames, 'stable');
            
            % don't pass the new options to the parent, or the
            % 'DefaultShape' option ( we set this to another value)
            nopass_list = [ { 'DefaultShape'; }; ...
                            C ];
            
        end
        
    end
    
end