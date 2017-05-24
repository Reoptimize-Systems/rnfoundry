classdef inline < mbdyn.pre.twoNodeJoint
    % Forces a point relative to the second node to move along a line
    % attached to the first node
    %
    % A point, optionally offset by relative_offset from the position of
    % node node_2_label, slides along a line that passes through a point
    % that is rigidly offset by relative_line_position from the position of
    % node_1_label, and is directed as direction 3 of relative orientation.
    
    
    properties (GetAccess = public, SetAccess = protected)
        
        relativeLinePosition;
        relativeOrientation;
        relativeOffset;
        linePositionReference;
        orientationReference;
        offsetReference;
        
    end
    
    methods
        
        function self = inline (node1, node2, varargin)
            
            options.RelativeLinePosition =  [];
            options.RelativeOrientation =  [];
            options.RelativeOffset = [];
            options.LinePositionReference = 'node';
            options.OrientationReference = 'node';
            options.OffsetReference = 'node';
            
            options = parse_pv_pairs (options, varargin);
            
            % call the superclass constructor
            self = self@mbdyn.pre.twoNodeJoint (node1, node2);
            
            self.type = 'in line';
            
            if ~isempty (options.RelativeLinePosition)
                self.relativeLinePosition = self.checkJointPositionOffset ({options.LinePositionReference, options.RelativeLinePosition});
                self.linePositionReference = options.LinePositionReference;
            else
                self.relativeLinePosition = options.RelativeLinePosition;
            end
            
            if ~isempty (options.RelativeOrientation)
                self.relativeOrientation = self.checkJointOrientationOffset ({options.OrientationReference, options.RelativeOrientation});
                self.orientationReference = options.OrientationReference;
            else
                self.relativeOrientation = options.RelativeOrientation;
            end
            
            if ~isempty (options.RelativeOffset)
                self.relativeOffset = self.checkJointPositionOffset ({options.OffsetReference, options.RelativeOffset});
                self.offsetReference = options.OffsetReference;
            else
                self.relativeOffset = options.RelativeOffset;
            end
            
        end
        
        function str = generateOutputString (self)
            
            str = generateOutputString@mbdyn.pre.twoNodeJoint(self);
            
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