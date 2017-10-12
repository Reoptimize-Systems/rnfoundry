classdef inPlane < mbdyn.pre.twoNodeJoint
    % This joint forces a point relative to the second node to move in a
    % plane attached to the first node
    %
    % A point, optionally offset by relative_offset from the position of
    % node 2, slides on a plane that passes through a point that is rigidly
    % offset by relative_plane_position from the position of node 1, and is
    % perpendicular to relative_normal. The vector relative_normal is
    % internally normalized to unity by mbdyn.
    %
    
    properties (GetAccess = public, SetAccess = protected)
        
        relativeNormal;
        relativePlanePosition;
        relativePlanePositionReference;
        relativeOffset;
        relativeOffsetReference;
        
    end
    
    methods
        
        function self = inPlane (node1, node2, relative_normal, varargin)
            
            options.RelativePlanePosition =  [];
            options.RelativePlanePositionReference = 'node';
            options.RelativeOffset = [];
            options.RelativeOffsetReference = 'node';
            
            options = parse_pv_pairs (options, varargin);
            
            % call the superclass constructor
            self = self@mbdyn.pre.twoNodeJoint (node1, node2);
            
            self.type = 'in plane';
            
            self.checkCartesianVector (relative_normal);
            
            self.relativeNormal = relative_normal;
            
            if ~isempty (options.RelativePlanePosition)
                self.checkCartesianVector (options.RelativePlanePosition);
                self.relativePlanePosition = self.checkJointPositionOffset ({options.RelativePlanePositionReference, options.RelativePlanePosition});
                self.relativePlanePositionReference = options.RelativePlanePositionReference;
            else
                self.relativePlanePosition = [];
            end
            
            if ~isempty (options.RelativeOffset)
                self.checkCartesianVector (options.RelativeOffset);
                self.relativeOffset = self.checkJointPositionOffset ({options.RelativeOffsetReference, options.RelativeOffset});
                self.relativeOffsetReference = options.RelativeOffsetReference;
            else
                self.relativeOffset = [];
            end
            
        end
        
        function str = generateOutputString (self)
            
            str = generateOutputString@mbdyn.pre.twoNodeJoint(self);
            
            str = self.addOutputLine (str, sprintf('%d', self.node1.label), 2, true, 'node 1 label');
            
            if isempty (self.relativePlanePosition)
                out = {'null'};
            else
                out = self.makeCellIfNot (self.relativePlanePosition);
            end
            str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 3, true);
            
            str = self.addOutputLine (str, self.commaSepList (self.relativeNormal), 3, true);
            
            addcomma = ~isempty (self.relativeOffset);
            str = self.addOutputLine (str, sprintf('%d', self.node2.label), 2, addcomma, 'node 2 label');
            
            if ~isempty (self.relativeOffset)
                out = self.makeCellIfNot (self.relativeOffset);
                str = self.addOutputLine (str, self.commaSepList ('offset', out{:}), 3, false);
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