classdef sphericalPin < mbdyn.pre.singleNodeOffsetJoint
    
    
    properties (GetAccess = public, SetAccess = protected)
 
        absolutePinPosition;
        absolutePinPositionReference;
        absolutePinOrientation;
        absolutePinOrientationReference;
        
    end
    
    methods
        
        function self = sphericalPin (node, absolute_pin_position, relative_offset, varargin)
            % sphericalPin constructor
            %
            % Syntax
            %
            %  rh = sphericalPin (node1)
            %  rh = sphericalPin (..., 'Parameter', value)
            %
            % Input
            %
            %  node - mbdyn.pre.structuralNode (or derived class) object
            %    representing the node the joint connects
            %
            %  absolute_pin_position - 
            %
            %
            % Additional arguments can be supplied as parameter-value
            % pairs. Available options are:
            %
            %  'PositionReference' - by default the position provided in
            %    absolute_pin_position is relative to the global reference
            %    frame. An alternative reference frame can be provided
            %    using this argument. Possible value for this are:
            %      'node'          : the default behaviour
            %      'global'        : the global reference frame
            %      'local'         : same as 'node'  
            %
            %  'AbsoluteOrientation' - mbdyn.pre.orientmat object
            %    containing the orientation of the joint relative to the
            %    gloabl frame. To provide an alternative reference you can
            %    use the optional Orientation1Reference parameter (see
            %    below)
            %
            %  'AbsoluteOrientationReference' - string containing a 
            %    reference for the orientation in AbsoluteOrientation, can
            %    be one of 'node', 'local' (equivalent to 'node') and
            %    'global'. Defaut is 'node'. 
            %  
            % Output
            %
            %  rh - mbdyn.pre.sphericalPin object
            %
            %
            
            [ options, nopass_list ] = mbdyn.pre.sphericalPin.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs ( options, nopass_list);
            
            % call the superclass constructor
            self = self@mbdyn.pre.singleNodeOffsetJoint (node, ...
                        'RelativeOffset', relative_offset, ...
                        pvpairs{:} );
                    
            allowedposrefstrs = {'global', 'node', 'local'};
            allowedorientrefstrs = {'global', 'node', 'local'};
            self.checkAllowedStringInputs ( options.PositionReference, allowedposrefstrs, true, 'PositionReference');
            self.checkAllowedStringInputs ( options.AbsoluteOrientationReference, allowedorientrefstrs, true, 'AbsoluteOrientationReference');
            self.checkOrientationMatrix (options.AbsoluteOrientation, true, 'AbsoluteOrientation');
            self.checkCartesianVector (absolute_pin_position, true, 'absolute_pin_position');
            
            self.absolutePinPosition = absolute_pin_position;
            self.absolutePinPositionReference = options.PositionReference;
            self.absolutePinOrientation = options.AbsoluteOrientation;
            self.absolutePinOrientationReference = options.AbsoluteOrientationReference;
            
            self.type = 'spherical pin';
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for sphericalPin joint
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (sp)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  sp - mbdyn.pre.sphericalPin object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %

            str = generateMBDynInputString@mbdyn.pre.singleNodeOffsetJoint (self, true);
            
            addcomma = ~isempty (self.absolutePinOrientation);
                    
            str = self.addOutputLine (str, self.commaSepList ('position', 'reference', self.absolutePinPositionReference, self.absolutePinPosition), 2, addcomma);

            if ~isempty (self.absolutePinOrientation)
                str = self.addOutputLine (str, self.commaSepList ('orientation', 'reference', self.absolutePinOrientationReference, self.absolutePinOrientation), 2, false);
            end
            
            str = self.addOutputLine (str, ';', 1, false, sprintf ('end %s', self.type));
            
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
            
            node1pos = self.node.absolutePosition;
            jref = self.reference ();
            jpos = jref.pos ();
                
            if ~self.needsRedraw
                % always have to redraw the line connecting the two points.
                % This changes shape, so we can't just transform the line
                % object
                delete (self.shapeObjects{1})
                self.shapeObjects{1} =  line ( self.drawAxesH, ...
                                               [ node1pos(1), jpos(1) ], ...
                                               [ node1pos(2), jpos(2) ], ...
                                               [ node1pos(3), jpos(3) ], ...
                                               'Color', self.drawColour );
                                       
            end
            
            if isempty (self.shapeObjects) ...
                    || self.needsRedraw
                % a full redraw is needed (and not just a modification of
                % transform matrices for the objects).
                
                % delete the current patch object
                self.deleteAllDrawnObjects ();
                
                self.shapeObjects = { line( self.drawAxesH, ...
                                            [ node1pos(1), jpos(1) ], ...
                                            [ node1pos(2), jpos(2) ], ...
                                            [ node1pos(3), jpos(3) ], ...
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
%             
%             om = self.absoluteJointOrientation;
%             
%             M = [ om.orientationMatrix, self.absoluteJointPosition; ...
%                   0, 0, 0, 1 ];
%                   
%             set ( self.transformObject, 'Matrix', M );
%             
        end
        
    end
    
    methods (Static)
        
        function [ options, nopass_list ] = defaultConstructorOptions ()
            
            options = mbdyn.pre.singleNodeOffsetJoint.defaultConstructorOptions ();
            
            parentfnames = fieldnames (options);
            
            options.AbsoluteOrientation = [];
            options.AbsoluteOrientationReference = 'node';
            options.PositionReference = 'global';
        
            allfnames = fieldnames (options);
            
            nopass_list = [ {'RelativeOffset'}; setdiff(allfnames, parentfnames)];
            
        end
        
    end
    
end