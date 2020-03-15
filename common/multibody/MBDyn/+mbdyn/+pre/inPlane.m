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
        relativeNormalReference;
        relativePlanePosition;
        relativePlanePositionReference;
        relativeOffset;
        relativeOffsetReference;
        
    end
    
    methods
        
        function self = inPlane (node1, node2, relative_normal, varargin)
        % mbdyn.pre.inPlane constructor
        %
        % Syntax
        %
        % ipo = mbdyn.pre.inPlane (node1, node2, relative_normal)
        % ipo = mbdyn.pre.inPlane (..., 'Parameter', Value)
        %
        % Description
        %
        % mbdyn.pre.inPlane forces a point relative to the second node to
        % move in a plane attached to the first node.
        %
        % A point, optionally offset by relative_offset from the position
        % of node 2, slides on a plane that passes through a point that is
        % rigidly offset by relative_plane_position from the position of
        % node 1, and is perpendicular to relative_normal. The vector
        % relative_normal defining the plane orientation is internally
        % normalized to unity by mbdyn.
        %
        % Input
        %
        %  node1 - node to which the plane will be attached
        %
        %  node2 - node which will be constrained to move on the plane.
        %
        %  relative_normal - optional (3 x 1) vector giving the
        %   normal of the plane specified (by default) in the reference
        %   frame of node 1. An alternative reference frame can be
        %   specified using the 'RelativeNormalReference' option described
        %   below.
        %
        % Addtional arguments may be supplied as parameter-value pairs. The
        % available options are:
        %
        %  'RelativeNormalReference' - optional string containing an
        %    alternative reference for the relative_normal. Can be one
        %    of 'global', 'node', 'local', other node', 'other
        %    position'.
        %
        %  'RelativePlanePosition' - optional (3 x 1) vector giving the
        %    relative offset of the plane from node 1.
        %
        %  'RelativePlanePositionReference' - optional character vector 
        %    containing an alternative reference for the
        %    RelativePlanePosition. Can be one of 'global', 'node',
        %    'local', other node', 'other position'.
        %
        %  'RelativeOffset' - optional (3 x 1) vector giving the
        %    relative offset from node 2 of the point constrained to move
        %    on the plane.
        %
        %  'RelativeOffsetReference' - optional string containing an
        %    alternative reference for the RelativeOffset. Can be one
        %    of 'global', 'node', 'local', other node', 'other
        %    position'.
        %
        % Output
        %
        %  ipo - mbdyn.pre.inPlane object
        %
        % See Also: mbdyn.pre.inLine
        %

            [ options, nopass_list ] = mbdyn.pre.inPlane.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs ( options, nopass_list);
            
            % call the superclass constructor
            self = self@mbdyn.pre.twoNodeJoint (node1, node2, pvpairs{:});
            
            self.type = 'in plane';
            
            self.checkCartesianVector (relative_normal);
            
            self.relativeNormal = self.checkJointPositionOffset ({options.RelativeNormalReference, relative_normal});
            self.relativeNormalReference = options.RelativeNormalReference;
            
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
        
        function str = generateMBDynInputString (self)
            
            str = generateMBDynInputString@mbdyn.pre.twoNodeJoint(self);
            
            str = self.addOutputLine (str, sprintf('%d', self.node1.label), 2, true, self.nodeLabelComment (self.node1));
            
            if isempty (self.relativePlanePosition)
                out = {'null'};
            else
                out = self.makeCellIfNot (self.relativePlanePosition);
            end
            str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 3, true);
            
            out = self.makeCellIfNot (self.relativeNormal);
            str = self.addOutputLine (str, self.commaSepList (out{:}), 3, true);
            
            addcomma = ~isempty (self.relativeOffset);
            str = self.addOutputLine (str, sprintf('%d', self.node2.label), 2, addcomma, self.nodeLabelComment (self.node2));
            
            if ~isempty (self.relativeOffset)
                out = self.makeCellIfNot (self.relativeOffset);
                str = self.addOutputLine (str, self.commaSepList ('offset', out{:}), 3, false);
            end
            
            str = self.addOutputLine (str, ';', 1, false, sprintf('end %s', self.type));
            
            self.addRegularization (str);
            
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
    
    methods (Static)
        
        function [ options, nopass_list ] = defaultConstructorOptions ()
            
            options = mbdyn.pre.twoNodeJoint.defaultConstructorOptions ();
            
            parentfnames = fieldnames (options);
            
            options.RelativePlanePosition =  [];
            options.RelativePlanePositionReference = 'node';
            options.RelativeOffset = [];
            options.RelativeOffsetReference = 'node';
            options.RelativeNormalReference = 'node';
            
            allfnames = fieldnames (options);
            
            C = setdiff (allfnames, parentfnames);
            
            nopass_list = C;
            
        end
        
    end
    
end