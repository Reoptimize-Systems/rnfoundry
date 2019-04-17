classdef sphericalHinge < mbdyn.pre.twoNodeOffsetJoint
    
    
    properties (GetAccess = public, SetAccess = protected)
%         frictionRadius;
%         frictionModel;
%         frictionShapeFcn;
%         preload;
    end
    
    methods
        
        function self = sphericalHinge (node1, node2, varargin)
            % sphericalHinge constructor
            %
            % Syntax
            %
            %  rh = sphericalHinge (node1, node2)
            %  rh = sphericalHinge (..., 'Parameter', value)
            %
            % Input
            %
            %  node1 - mbdyn.pre.structuralNode (or derived class) object
            %    representing to first node the joint connects
            %
            %  node2 - mbdyn.pre.structuralNode (or derived class) object
            %    representing to second node the joint connects
            %
            % Additional arguments can be supplied as parameter-value
            % pairs. Available options are:
            %
            %  'RelativeOffset1' - (3 x 1) vector containing the offset of 
            %    the ball joint relative to the first node. To provide an
            %    alternative reference you can use the optional
            %    Offset1Reference parameter (see below). Default is null
            %    (no offset) if not supplied.
            %
            %  'RelativeOffset2' - (3 x 1) vector containing the offset of 
            %    the ball joint relative to the second node. To provide an
            %    alternative reference you can use the optional
            %    Offset1Reference parameter (see below). Default is null
            %    (no offset) if not supplied.
            %
            %  'Offset1Reference' - by default the positions provided in
            %    position1 and position2 are relaive to the respective
            %    nodes in their reference frame. An alternative reference
            %    frame can be provided using this argument. Possible
            %    value for this are: 
            %      'node'          : the default behaviour
            %      'global'        : the global reference frame
            %      'other node'    : the frame of the other node the joint  
            %                        is attached to
            %      'other position': a relative position in the other 
            %                        node's reference frame, with respect 
            %                        to the relative position already 
            %                        specified for the other node
            %
            %  'Offset2Reference' - same as Offset1Reference, but for the
            %    second node
            %
            %  'RelativeOrientation1' - mbdyn.pre.orientmat object
            %    containing the orientation of the joint relative to the
            %    first node. To provide an alternative reference you can
            %    use the optional Orientation1Reference parameter (see
            %    below)
            %
            %  'RelativeOrientation2' - mbdyn.pre.orientmat object
            %    containing the orientation of the joint relative to the
            %    second node. To provide an alternative reference you can
            %    use the optional Orientation2Reference parameter (see
            %    below)
            %
            %  'Orientation1Reference' - string containing a reference for
            %    the orientation in RelativeOrientation1, can be one of
            %    'node', 'local' (equivalent to 'node'), 'other node',
            %    'other orientation' and 'global'. Defaut is 'node'. See
            %    Offset1Reference above for more information.
            %
            %  'Orientation2Reference' - string containing a reference for
            %    the orientation in RelativeOrientation2, can be one of
            %    'node', 'local' (equivalent to 'node'), 'other node',
            %    'other orientation' and 'global'. Defaut is 'node'. See
            %    Offset1Reference above for more information.
            %
            %  
            % Output
            %
            %  rh - mbdyn.pre.sphericalHinge object
            %
            %
            
            [ options, nopass_list ] = mbdyn.pre.sphericalHinge.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs ( options, nopass_list);
            
            % call the superclass constructor
            self = self@mbdyn.pre.twoNodeOffsetJoint (node1, node2, ...
                        pvpairs{:} );
            
            self.type = 'spherical hinge';
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for sphericalHinge joint
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (rh)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  rh - mbdyn.pre.sphericalHinge object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %

            str = generateMBDynInputString@mbdyn.pre.twoNodeOffsetJoint (self, false);
            
            str = self.addOutputLine (str, ';', 1, false, sprintf('end %s', self.type));
            
        end
        
        function draw (self, varargin)
            
            options.AxesHandle = [];
            options.ForceRedraw = false;
            options.Mode = 'solid';
            
            options = parse_pv_pairs (options, varargin);
            
%             draw@mbdyn.pre.twoNodeOffsetJoint ( self, ...
%                 'AxesHandle', options.AxesHandle, ...
%                 'ForceRedraw', options.ForceRedraw, ...
%                 'Mode', options.Mode );

            if options.ForceRedraw
                self.needsRedraw = true;
            end
            
            self.checkAxes (options.AxesHandle);
            
            node1pos = self.node1.absolutePosition;
            jref = self.reference ();
            jpos = jref.pos ();
            node2pos = self.node2.absolutePosition;
                
            if ~self.needsRedraw
                % always have to redraw the line connecting the two points.
                % This changes shape, so we can't just transform the line
                % object
                delete (self.shapeObjects{1})
                self.shapeObjects{1} =  line ( self.drawAxesH, ...
                                               [ node1pos(1), jpos(1), node2pos(1) ], ...
                                               [ node1pos(2), jpos(2), node2pos(2) ], ...
                                               [ node1pos(3), jpos(3), node2pos(3) ], ...
                                               'Color', self.drawColour );
                                       
            end
            
            if isempty (self.shapeObjects) ...
                    || self.needsRedraw
                % a full redraw is needed (and not just a modification of
                % transform matrices for the objects).
                
                % delete the current patch object
                self.deleteAllDrawnObjects ();
                
                self.shapeObjects = { line( self.drawAxesH, ...
                                            [ node1pos(1), jpos(1), node2pos(1) ], ...
                                            [ node1pos(2), jpos(2), node2pos(2) ], ...
                                            [ node1pos(3), jpos(3), node2pos(3) ], ...
                                            'Color', self.drawColour ), ...
                                      line( self.drawAxesH, ...
                                            [ 0, 0 ], ...
                                            [ 0, 0 ], ...
                                            [ -self.shapeParameters(1)/2, self.shapeParameters(3)/2 ], ...
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

%             self.setTransform ();
            
        end
        
    end
    
    methods (Access = protected)
        
        function setTransform (self)
            
%             om = self.absoluteJointOrientation;
%             
%             M = [ om.orientationMatrix, self.absoluteJointPosition; ...
%                   0, 0, 0, 1 ];
%                   
%             set ( self.transformObject, 'Matrix', M );
            
        end
        
    end
    
    methods (Static)
        
        function [ options, nopass_list ] = defaultConstructorOptions ()
            
            options = mbdyn.pre.twoNodeOffsetJoint.defaultConstructorOptions ();
            
            nopass_list = {};
            
        end
        
    end
    
end