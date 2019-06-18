classdef inline < mbdyn.pre.twoNodeJoint

    properties (GetAccess = public, SetAccess = protected)
        
        relativeLinePosition;
        relativeOrientation;
        relativeOffset;
        linePositionReference;
        orientationReference;
        offsetReference;
        
    end
    
    properties (GetAccess = public, SetAccess = protected)
        
        lineTransformObj;
        
    end
    
    methods
        
        function self = inline (node1, node2, varargin)
            % constructor for inline joint
            %
            %
            % Syntax
            %
            % inlobj = inline (node1, node2)
            % inlobj = inline (..., 'Parameter', value)
            %
            % Description
            %
            % Joint which forces a point relative to the second node to
            % move along a line attached to the first node.
            %
            % A point, optionally offset from the position of the second
            % node, slides along a line that passes through a point that is
            % rigidly offset from the position of node 1, and is directed
            % as direction 3 of an orientation relative to node 1. See the
            % description of the joint in the MBDyn manual for more
            % information.
            %
            % Input
            %
            %  node1 - mbdyn.pre.structuralNode6Dof object
            %
            %  node2 - mbdyn.pre.structuralNode6Dof object
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'RelativeLinePosition' - optional (3 x 1) vector giving the
            %    relative offset of the line from node 1.
            %
            %  'RelativeOrientation' - optional mbdyn.pre.orientmat object
            %    giving the relative orientation of the line in the frame
            %    of node 1. The line is parallel to axis three of the
            %    resulting orientation.
            %
            %  'RelativeOffset' - optional (3 x 1) vector giving the
            %    relative offset from node 2 of the point constrained to
            %    move along the line.
            %
            %  'LinePositionReference' - optional string containing an
            %    alternative reference for the RelativeLinePosition. Can be
            %    one of 'global', 'node', 'local', other node', 'other
            %    position'.
            %
            %  'OrientationReference' - optional string containing an
            %    alternative reference for the RelativeOrientation. Can be
            %    one of 'global', 'node', 'local', other node', 'other
            %    orientation'.
            %
            %  'OffsetReference' - optional string containing an
            %    alternative reference for the RelativeOffset. Can be one
            %    of 'global', 'node', 'local', other node', 'other
            %    position'.
            %
            % Output
            %
            %  inlobj - mbdyn.pre.inline object
            %
            %
            %
            % See Also: 
            %
            
            [options, nopass_list] = mbdyn.pre.inline.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            
            % call the superclass constructor
            self = self@mbdyn.pre.twoNodeJoint (node1, node2, pvpairs{:}, 'DefaultShape', 'none');
            
            self.type = 'in line';
            
            if ~isempty (options.RelativeLinePosition)
                self.relativeLinePosition = self.checkJointPositionOffset ({options.LinePositionReference, options.RelativeLinePosition});
            else
                self.relativeLinePosition = options.RelativeLinePosition;
            end
            
            if ~isempty (options.RelativeOrientation)
                self.relativeOrientation = self.checkJointOrientationOffset ({options.OrientationReference, options.RelativeOrientation});
            else
                self.relativeOrientation = options.RelativeOrientation;
            end
            
            if ~isempty (options.RelativeOffset)
                self.relativeOffset = self.checkJointPositionOffset ({options.OffsetReference, options.RelativeOffset});
            else
                self.relativeOffset = options.RelativeOffset;
            end
            
            self.linePositionReference = options.LinePositionReference;
            self.orientationReference = options.OrientationReference;
            self.offsetReference = options.OffsetReference;
            
            self.setSize ();
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = generateMBDynInputString@mbdyn.pre.twoNodeJoint(self);
            
            str = self.addOutputLine (str, sprintf('%d', self.node1.label), 2, true, 'node 1 label');
            
            if ~isempty (self.relativeLinePosition)
                out = self.makeCellIfNot (self.relativeLinePosition);
                str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 3, true, 'relative line position');
                
                if ~isempty (self.relativeOrientation)
                    out = self.makeCellIfNot (self.relativeOrientation);
                    str = self.addOutputLine (str, self.commaSepList ('orientation', out{:}), 3, true, 'relative orientation');
                end
            
            end
            
            
            
            addcomma = ~isempty (self.relativeOffset);
            str = self.addOutputLine (str, sprintf('%d', self.node2.label), 2, addcomma, 'node 2 label');
            
            if ~isempty (self.relativeOffset)
                out = self.makeCellIfNot (self.relativeOffset);
                str = self.addOutputLine (str, self.commaSepList ('offset', out{:}), 3, false);
            end
            
            str = self.addOutputLine (str, ';', 1, false, sprintf('end %s', self.type));
            
        end
        
        
        function setSize (self, varargin)
            % set the size of the element in plots
            %
            % Syntax
            %
            % setSize (el)
            % setSize (..., 'Parameter', value)
            %
            % Description
            %
            % setSize is used to set the sizes of the inline joint shapes
            % for plotting the element in a figure.
            %
            % Input
            %
            %  el - mbdyn.pre.element object
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'LineLength' - length of the line the point slides on in the
            %    plot
            %
            %  'CylinderRadius' - radius of the cylinder representing the
            %    point on the line
            %
            %  'CylinderLength' - axial length of the cylinder representing
            %    the point on the line
            %
            % See Also: 
            %

            options.LineLength = 1;
            options.CylinderRadius = [];
            options.CylinderLength = [];

            options = parse_pv_pairs (options, varargin);

            self.checkNumericScalar (options.LineLength, true, 'LineLength');

            if isempty (options.CylinderRadius)
                options.CylinderRadius = 0.05 * options.LineLength;
            end
            if isempty (options.CylinderLength)
                options.CylinderLength = 0.1 * options.LineLength;
            end

            self.checkNumericScalar (options.CylinderRadius, true, 'CylinderRadius');
            self.checkNumericScalar (options.CylinderLength, true, 'CylinderLength');

            assert (options.LineLength > 0, 'LineLength must be greater than zero');
            assert (options.CylinderRadius > 0, 'CylinderRadius must be greater than zero');
            assert (options.CylinderLength > 0, 'CylinderLength must be greater than zero');

            self.shapeParameters(1) = options.LineLength;
            self.shapeParameters(2) = options.CylinderRadius;
            self.shapeParameters(3) = options.CylinderLength;

            % set the shapedata to empty so it is recreated with the new
            % sizes when draw is next called
            self.shapeData = [];
            self.needsRedraw = true;
            
        end
        
        
        function draw (self, varargin)
            
            options.AxesHandle = [];
            options.ForceRedraw = false;
            options.Mode = [];
            
            options = parse_pv_pairs (options, varargin);
              
            if options.ForceRedraw
                self.needsRedraw = true;
            end
            
            self.checkAxes (options.AxesHandle);
            
            if isempty (self.shapeObjects) ...
                    || self.needsRedraw
                % a full redraw is needed (and not just a modification of
                % transform matrices for the objects).
                
                if isempty (self.lineTransformObj) || ~ishghandle (self.lineTransformObj)
                    self.lineTransformObj = hgtransform (self.drawAxesH);
                end
                
                self.shapeData = self.makeAnnularCylinderShape ( self.shapeParameters(2), ...
                                                                 self.shapeParameters(2)/2, ...
                                                                 self.shapeParameters(3) );
                
                
                
                % delete the objects
                self.deleteAllDrawnObjects ();
                
                % make the line
                self.shapeObjects = { line( self.drawAxesH, ...
                                            [ 0, 0 ], ...
                                            [ 0, 0 ], ...
                                            [ -self.shapeParameters(1)/2, self.shapeParameters(1)/2 ], ...
                                            'Color', self.drawColour, ...
                                            'LineStyle', '--', ...
                                            'Parent', self.lineTransformObj ), ...
                                    };
                                

                
                % the point 
                for ind = 1:numel (self.shapeData)
                    
                    self.shapeData{ind}.FaceColor = self.drawColour;
                    self.shapeData{ind}.FaceAlpha = 0.5;
                    self.shapeData{ind}.EdgeColor = self.drawColour;
                    self.shapeData{ind}.FaceLighting = 'Gouraud';
                    self.shapeData{ind}.AmbientStrength = 0.15;
                    self.shapeData{ind}.Parent = self.transformObject;
                    
                    self.shapeObjects = [ self.shapeObjects, ...
                                          { patch( self.drawAxesH, ...
                                                   self.shapeData{ind} ) } ...
                                        ];
                end

                                     
                self.needsRedraw = false;

                
            end
            
            self.setTransform ();


        end
        
        function abspos = lineAbsolutePosition (self)
            % gets the position of the clamp in the global frame

            abspos = offset2AbsolutePosition ( self, ...
                                               self.relativeLinePosition{3}, ...
                                               self.linePositionReference, ...
                                               1 );
            
        end
        
        function absorientm = lineAbsoluteOrientation (self)
            % gets the orientation of the clamp in the global frame
            
            absorientm = orient2AbsoluteOrientation ( ...
                                               self, ...
                                               self.relativeOrientation{3}, ...
                                               self.orientationReference, ...
                                               1 );
            
        end
        
        function abspos = pointAbsolutePosition (self)
            % gets the position of the sliding point in the global frame
            
            if isempty (self.relativeOffset)
                offset = [0;0;0];
            else
                offset = self.relativeOffset;
            end
            
            abspos = offset2AbsolutePosition ( self, ...
                                               offset , ...
                                               self.offsetReference, ...
                                               2 );
            
        end
        
        
    end
    
    methods (Access = protected)
        
        function setTransform (self)
            
            linepos = lineAbsolutePosition (self);

            lineorientm = lineAbsoluteOrientation (self);

            pointpos = pointAbsolutePosition (self);
                     
%             % modify the line data
%             set ( self.shapeObjects{1}, ...
%                   'XData', [ lineends(1,1), lineends(1,2) ], ...
%                   'YData', [ lineends(2,1), lineends(2,2) ], ...
%                   'ZData', [ lineends(3,1), lineends(3,2) ] );
              
            % set the transform object which controls the line
            % location and orientation
            M = [ lineorientm, linepos; ...
                  0, 0, 0, 1 ];
                  
            set ( self.lineTransformObj, 'Matrix', M );
            
            % set the transform object which controls the point cylinder
            % location and orientation
            M = [ lineorientm, pointpos; ...
                  0, 0, 0, 1 ];
                  
            set ( self.transformObject, 'Matrix', M );
            
            
        end
        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.twoNodeJoint.defaultConstructorOptions ();
            
            parentfnames = fieldnames (options);
            
            % add default options common to all inline objects
            options.RelativeLinePosition =  'null';
            options.RelativeOrientation = mbdyn.pre.orientmat ('eye');
            options.RelativeOffset = [];
            options.LinePositionReference = 'node';
            options.OrientationReference = 'node';
            options.OffsetReference = 'node';
            
            allfnames = fieldnames (options);
            
            nopass_list = [ setdiff(allfnames, parentfnames); {'DefaultShape'}];
            
        end
        
    end
    
end