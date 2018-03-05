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
            % compatible: use the orientation keyword to assign the joint
            % initial orientation.
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

            options.RelativeOrientation1 =  [];
            options.RelativeOrientation2 =  [];
            options.Orientation1Reference = 'node';
            options.Orientation2Reference = 'node';
            
            options = parse_pv_pairs (options, varargin);
            
            % call the superclass constructor
            self = self@mbdyn.pre.twoNodeJoint (node1, node2);
            
            self.type = 'prismatic';
            
            if ~isempty (self.relativeOrientation1)
                self.relativeOrientation1 = self.checkJointOrientationOffset ({options.Orientation1Reference, ptions.RelativeOrientation1});
                self.orientation1Reference = options.Orientation1Reference;
            else
                self.relativeOrientation1 = [];
            end
            
            if ~isempty (self.relativeOrientation2)
                self.relativeOrientation2 = self.checkJointOrientationOffset ({options.Orientation2Reference, options.RelativeOrientation2});
                self.orientation2Reference = options.Orientation2Reference;
            else
                self.relativeOrientation2 = [];
            end
            
        end
        
        function str = generateOutputString (self)
            
            str = generateOutputString@mbdyn.pre.twoNodeJoint(self);
            
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
            
%             options.AxesHandle = [];
%             options.ForceRedraw = false;
%             options.Mode = 'solid';
%             
%             options = parse_pv_pairs (options, varargin);
%             
%             draw@mbdyn.pre.element ( self, ...
%                 'AxesHandle', options.AxesHandle, ...
%                 'ForceRedraw', options.ForceRedraw, ...
%                 'Mode', options.Mode );
% 
%             self.setTransform ();
            
        end
        
    end
    
    methods (Access = protected)
        
        function setTransform (self)
            
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
            
        end
        
    end
    
end