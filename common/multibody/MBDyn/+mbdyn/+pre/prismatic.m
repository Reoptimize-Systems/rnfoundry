classdef prismatic < mbdyn.pre.twoNodeJoint
   
    
    properties (GetAccess = public, SetAccess = protected)
        
        relativeOrientation1;
        relativeOrientation2;
        
        orientation1Reference;
        orientation2Reference;
        
    end
    
    methods
        
        function self = prismatic (node1, node2, varargin)
            % prismatic joint constructor
            %
            % Syntax
            %
            % pjnt = prismatic (node1, node2)
            % pjnt = prismatic (.., 'Parameter', value)
            %
            % Description
            %
            % Constrains the relative orientation of two nodes, so that
            % their orientations remain parallel. The relative position is
            % not constrained. The initial orientation of the joint must be
            % compatible: use the RelativeOrientation1 option to assign the
            % joint initial orientation.
            %
            % Input
            %
            %  node1 - mbdyn.pre.structuralNode6dof object
            %
            %  node2 - mbdyn.pre.structuralNode6dof object
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'RelativeOrientation1' - 
            %
            %  'RelativeOrientation2' - 
            %
            %  'Orientation1Reference' - 
            %
            %  'Orientation2Reference' - 
            %
            % Output
            %
            %  pjnt - mbdyn.pre.prismatic object
            %
            %
            %
            % See Also: 
            %

            [options, nopass_list] = mbdyn.pre.prismatic.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            % call the superclass constructor
            self = self@mbdyn.pre.twoNodeJoint (node1, node2, pvpairs{:});
            
            self.type = 'prismatic';
            
            if ~isempty (options.RelativeOrientation1)
                self.relativeOrientation1 = self.checkJointOrientationOffset ({options.Orientation1Reference, options.RelativeOrientation1});
                self.orientation1Reference = options.Orientation1Reference;
            else
                self.relativeOrientation1 = [];
            end
            
            if ~isempty (options.RelativeOrientation2)
                self.relativeOrientation2 = self.checkJointOrientationOffset ({options.Orientation2Reference, options.RelativeOrientation2});
                self.orientation2Reference = options.Orientation2Reference;
            else
                self.relativeOrientation2 = [];
            end
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = generateMBDynInputString@mbdyn.pre.twoNodeJoint(self);
            
            str = self.addOutputLine (str, sprintf('%d', self.node1.label), 2, true, 'node 1 label');
            
            if ~isempty (self.relativeOrientation1)
                out = self.makeCellIfNot (self.relativeOrientation1);
                str = self.addOutputLine (str, self.commaSepList ('orientation', out{:}), 3, true);
            end
            
            addcomma = ~isempty (self.relativeOrientation2);
            str = self.addOutputLine (str, sprintf('%d', self.node2.label), 2, addcomma, 'node 2 label');
            
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
            
%             draw@mbdyn.pre.twoNodeOffsetJoint ( self, ...
%                 'AxesHandle', options.AxesHandle, ...
%                 'ForceRedraw', options.ForceRedraw, ...
%                 'Mode', options.Mode );

            if options.ForceRedraw
                self.needsRedraw = true;
            end
            
            self.checkAxes (options.AxesHandle);
            
            node1pos = self.node1.absolutePosition;
            node2pos = self.node2.absolutePosition;
                
            if ~self.needsRedraw
                % always have to redraw the line connecting the two points.
                % This changes shape, so we can't just transform the line
                % object
                delete (self.shapeObjects{1})
                self.shapeObjects{1} =  line ( self.drawAxesH, ...
                                               [ node1pos(1), node2pos(1) ], ...
                                               [ node1pos(2), node2pos(2) ], ...
                                               [ node1pos(3), node2pos(3) ], ...
                                               'Color', self.drawColour );
                                       
            end
            
            if isempty (self.shapeObjects) ...
                    || self.needsRedraw
                % a full redraw is needed (and not just a modification of
                % transform matrices for the objects).
                
                % delete the current objects
                self.deleteAllDrawnObjects ();
                radius = max ( 1.05*[self.node1.sx, self.node1.sy, self.node1.sz]/2 );
                
                % circle centre locations will be set by the node
                % transform matrices
                circ1_points = self.circlePoints3D ([0;0;0], [1,0,0], radius, 20);
                circ2_points = self.circlePoints3D ([0;0;0], [0,1,0], radius, 20);
                circ3_points = self.circlePoints3D ([0;0;0], [0,0,1], radius, 20);
                
                self.shapeObjects = { line( self.drawAxesH, ...
                                            [ node1pos(1), node2pos(1) ], ...
                                            [ node1pos(2), node2pos(2) ], ...
                                            [ node1pos(3), node2pos(3) ], ...
                                            'Color', self.drawColour ), ...
                                      line( self.drawAxesH, ...
                                            circ1_points(1,:), ...
                                            circ1_points(2,:), ...
                                            circ1_points(3,:), ...
                                            'Color', self.drawColour, ...
                                            'Parent', self.node1.transformObject ), ...
                                      line( self.drawAxesH, ...
                                            circ2_points(1,:), ...
                                            circ2_points(2,:), ...
                                            circ2_points(3,:), ...
                                            'Color', self.drawColour, ...
                                            'Parent', self.node1.transformObject ), ...
                                      line( self.drawAxesH, ...
                                            circ3_points(1,:), ...
                                            circ3_points(2,:), ...
                                            circ3_points(3,:), ...
                                            'Color', self.drawColour, ...
                                            'Parent', self.node1.transformObject ), ...
                                      line( self.drawAxesH, ...
                                            circ1_points(1,:), ...
                                            circ1_points(2,:), ...
                                            circ1_points(3,:), ...
                                            'Color', self.drawColour, ...
                                            'Parent', self.node2.transformObject ), ...
                                      line( self.drawAxesH, ...
                                            circ1_points(1,:), ...
                                            circ1_points(2,:), ...
                                            circ1_points(3,:), ...
                                            'Color', self.drawColour, ...
                                            'Parent', self.node2.transformObject ), ...
                                      line( self.drawAxesH, ...
                                            circ1_points(1,:), ...
                                            circ1_points(2,:), ...
                                            circ1_points(3,:), ...
                                            'Color', self.drawColour, ...
                                            'Parent', self.node2.transformObject ) ...
                                    };
                
                self.needsRedraw = false;
                
            end
            
            % normally we would set the transform matrix around here, but
            % we don't need to for the prismatic joint because the nodes
            % set the transform matrix and we are using the same transform
%             self.setTransform ();
            
        end
        
    end
    
    methods (Access = protected)
        
%         function setTransform (self)
%             
%             switch self.orientation1Reference
%                 
%                 case 'node'
%                     ref_orient_base = mbdyn.pre.reference (self.node1.absolutePosition, ...
%                                                            self.node1.absoluteOrientation, ...
%                                                            self.node1.absoluteVelocity, ...
%                                                            self.node1.absoluteAngularVelocity);
%                 case 'global'
%                     ref_orient_base = mbdyn.pre.globalref;
%                     
%             end
% 
% %                                         
%             ref_joint = mbdyn.pre.reference (self.relativeOffset1{end}, ...
%                             mbdyn.pre.orientmat ('orientation', self.relativeOrientation1{end}), ...
%                             [], ...
%                             [], ...
%                             'PositionParent', ref_pos_base, ...
%                             'OrientParent', ref_orient_base);
%             
%             M = [ ref_joint.orientm.orientationMatrix , ref_joint.pos; ...
%                   0, 0, 0, 1 ];
%             
%             % matlab uses different convention to mbdyn for rotation
%             % matrix
%             M = self.mbdynOrient2Matlab (M);
%                   
%             set ( self.transformObject, 'Matrix', M );
%             
%         end
        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.twoNodeJoint.defaultConstructorOptions ();
            
            parentfnames = fieldnames (options);
            
            % add default options common to all twoNodeOffsetJoint objects
            options.RelativeOrientation1 =  [];
            options.RelativeOrientation2 =  [];
            options.Orientation1Reference = 'node';
            options.Orientation2Reference = 'node';
            
            allfnames = fieldnames (options);
            
            nopass_list = setdiff (allfnames, parentfnames);
            
        end
        
    end
    
end