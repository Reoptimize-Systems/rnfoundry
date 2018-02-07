classdef revoluteHinge < mbdyn.pre.twoNodeOffsetJoint
    
    
    properties (GetAccess = public, SetAccess = protected)
        frictionRadius;
        frictionModel;
        preload;
    end
    
    methods
        
        function self = revoluteHinge (node1, node2, position1, position2, varargin)
            % revoluteHinge constructor
            %
            % Syntax
            %
            %  rh = revoluteHinge (node1, node2, position1, position2, 'Parameter', value)
            %
            % Input
            %
            %  node1 - mbdyn.pre.structuralNode (or derived class) object
            %    representing to first node the joint connects
            %
            %  node2 - mbdyn.pre.structuralNode (or derived class) object
            %    representing to second node the joint connects
            %
            %  position1 - (3 x 1) vector containing the offset of the
            %    joint relative to the first node. To provide an
            %    alternative reference you can use the optional
            %    Offset1Reference parameter (see below)
            %
            %  position2 - (3 x 1) vector containing the offset of the
            %    joint relative to the second node. To provide an
            %    alternative reference you can use the optional
            %    Offset1Reference parameter (see below)
            %
            % Additional arguments can be supplied as parameter-value
            % pairs. Available options are:
            %
            %  'Offset1Reference' - by default the positions provided in
            %    position1 and position2 are relaive to the respective
            %    nodes in their reference frame. An alternative reference
            %    frame can be provided using this argument. Possible
            %    value for this are: 
            %      'node'          : the default behaviour
            %      'global' -      : the global reference frame
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
            %  'RelativeOrientation1' - 
            %
            %  'RelativeOrientation2' - 
            %
            %  'Orientation1Reference' = 'node';
            %
            %  'Orientation2Reference' - 
            %
            %  'InitialTheta' - 
            %
            %  'Friction' - 
            %
            %  'Preload'' - 
            %
            %  'ShapeFunction' - 
            %  
            % Output
            %
            %  rh - mbdyn.pre.revoluteHinge object
            %
            %
            
            options.RelativeOrientation1 =  [];
            options.RelativeOrientation2 =  [];
            options.InitialTheta = [];
            options.FrictionRadius = [];
            options.FrictionModel = [];
            options.Preload = [];
            options.ShapeFunction = [];
            options.Offset1Reference = 'node';
            options.Offset2Reference = 'node';
            options.Orientation1Reference = 'node';
            options.Orientation2Reference = 'node';
            
            options = parse_pv_pairs (options, varargin);
            
            % call the superclass constructor
            self = self@mbdyn.pre.twoNodeOffsetJoint (node1, node2, ...
                        'RelativeOffset1', position1, ...
                        'RelativeOffset2', position2, ...
                        'RelativeOrientation1', options.RelativeOrientation1, ...
                        'RelativeOrientation2', options.RelativeOrientation2, ...
                        'Offset1Reference', options.Offset1Reference, ...
                        'Offset2Reference', options.Offset2Reference, ...
                        'Orientation1Reference', options.Orientation1Reference, ...
                        'Orientation2Reference', options.Orientation2Reference );
            
            if ~isempty (options.FrictionRadius)
                if isempty (options.FrictionModel)
                    error ('If supplying a friction radius, you must also supply a friction model');
                end
                self.checkNumericScalar (options.FrictionRadius, true, 'FrictionRadius')
                
                if ~isempty (options.Preload)
                    self.checkNumericScalar (options.Preload, true, 'Preload');
                    self.preload = options.Preload;
                end
            end
            
            if ~isempty (options.FrictionModel)
                if isempty (options.FrictionRadius)
                    error ('If supplying a friction model, you must also supply a friction radius');
                end
                assert (isa (options.FrictionModel, 'mbdyn.pre.frictionModel'), ...
                    'Supplied FrictionModel is not an mbdyn.pre.frictionModel object (or derived class)');
            end
            
            self.frictionRadius = options.FrictionRadius;
            self.frictionModel = options.FrictionModel;
            
            self.type = 'revolute hinge';
            
        end
        
        function str = generateOutputString (self)
            
            str = generateOutputString@mbdyn.pre.twoNodeJoint(self);
            
            str = self.addOutputLine (str, sprintf('%d', self.node1.label), 2, true, 'node 1 label');
            
%             out = self.makeCellIfNot (self.relativeOffset1);
            str = self.addOutputLine (str, self.commaSepList ('position', 'reference', self.offset1Reference, self.relativeOffset1), 3, true);
            
            if ~isempty (self.relativeOrientation1)
%                 out = self.makeCellIfNot (self.relativeOrientation1);
                str = self.addOutputLine (str, self.commaSepList ('orientation', 'reference', self.orientation1Reference, self.relativeOrientation1), 3, true);
            end
            
            str = self.addOutputLine (str, sprintf('%d', self.node2.label), 2, true, 'node 2 label');
            
%             out = self.makeCellIfNot (self.relativeOffset2);
            addcomma = ~isempty (self.relativeOrientation2);
            str = self.addOutputLine (str, self.commaSepList ('position', 'reference', self.offset2Reference, self.relativeOffset2), 3, addcomma);
            
            addcomma = ~isempty (self.frictionRadius);
            if ~isempty (self.relativeOrientation2)
%                 out = self.makeCellIfNot (self.relativeOrientation2);
                str = self.addOutputLine (str, self.commaSepList ('orientation', 'reference', self.orientation2Reference, self.relativeOrientation2), 3, addcomma);
            end
            
            if ~isempty (self.frictionRadius)
                str = self.addOutputLine (str, self.commaSepList ('friction', self.frictionRadius), 3, true);
                
                if ~isempty (self.preload)
                    str = self.addOutputLine (str, self.commaSepList ('preload', self.frictionRadius), 4, true);
                end
                
                str = self.addOutputLine (str, self.frictionModel.generateOutputString (), 4, false);
            end
            
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
                % always have to redraw line, can't just transform objects
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

%             self.setTransform ();
            
        end
        
    end
    
    methods (Access = protected)
        
%         function setTransform (self)
%             
%             om = self.absoluteJointOrientation;
%             
%             M = [ om.orientationMatrix, self.absoluteJointPosition; ...
%                   0, 0, 0, 1 ];
%                   
%             set ( self.transformObject, 'Matrix', M );
%             
%         end
        
    end
    
end