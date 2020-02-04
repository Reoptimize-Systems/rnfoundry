classdef revoluteRotation < mbdyn.pre.twoNodeJoint
    
    properties (GetAccess = public, SetAccess = protected)
        
        relativeOffset1;
        relativeOffset1Reference;
        relativeOrientation1;
        relativeOrientation1Reference;
            
        relativeOffset2;
        relativeOffset2Reference;
        relativeOrientation2;
        relativeOrientation2Reference;
        
    end
    
    
    methods
        
        function self = revoluteRotation (node1, node2, varargin)
            % constructor for revolute rotation joint
            %
            % Syntax
            %
            % rrjnt = revoluteRotation (node1, node2)
            % rrjnt = revoluteRotation (..., 'Parameter', value)
            %
            % Description
            %
            % revoluteRotation is a joint which allows the relative
            % rotation of two nodes about a given axis, which is axis 3 in
            % the reference systems defined by two orientations. The
            % relative position is not constrained.
            %
            % Input
            %
            %  node1 - mbdyn.pre.structuralNode object
            %
            %  node2 - mbdyn.pre.structuralNode object
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'RelativeOffset1' - 
            %
            %  'RelativeOffset1Reference' - 
            %
            %  'RelativeOrientation1' - 
            %
            %  'RelativeOrientation1Reference' - 
            %
            %  'RelativeOffset2' - 
            %
            %  'RelativeOffset2Reference' - 
            %
            %  'RelativeOrientation2' - 
            %
            %  'RelativeOrientation2Reference' - 
            %
            % Output
            %
            %  rrjnt - mbdyn.pre.revoluteRotation object
            %
            %
            %
            % See Also: 
            %

            options.RelativeOffset1 = [];
            options.RelativeOffset1Reference = 'node';
            options.RelativeOrientation1 =  [];
            options.RelativeOrientation1Reference = 'node';
            
            options.RelativeOffset2 =  [];
            options.RelativeOffset2Reference = 'node';
            options.RelativeOrientation2 =  [];
            options.RelativeOrientation2Reference = 'node';
            
            options = parse_pv_pairs (options, varargin);
            
            % call the superclass constructor
            self = self@mbdyn.pre.twoNodeJoint (node1, node2);
            
            self.type = 'revolute rotation';
            
            if ~isempty (options.RelativeOffset1)
                self.relativeOffset1 = self.checkJointPositionOffset ( {options.RelativeOffset1Reference, options.RelativeOffset1});
                self.relativeOffset1Reference = options.RelativeOffset1Reference;
            else
                self.relativeOffset1 = options.RelativeOffset1;
            end
            
            if ~isempty (options.RelativeOffset2)
                self.relativeOffset2 = self.checkJointPositionOffset ( {options.RelativeOffset2Reference, options.RelativeOffset2});
                self.relativeOffset2Reference = options.RelativeOffset2Reference;
            else
                self.relativeOffset2 = options.RelativeOffset2;
            end
            
            if ~isempty (options.RelativeOrientation1)
                self.relativeOrientation1 = self.checkJointOrientationOffset ( {options.RelativeOrientation1Reference, options.RelativeOrientation1});
                self.relativeOrientation1Reference = options.RelativeOrientation1Reference;
            else
                self.relativeOrientation1 = options.RelativeOrientation1;
            end
            
            if ~isempty (options.RelativeOrientation2)
                self.relativeOrientation2 = self.checkJointOrientationOffset ( {options.RelativeOrientation2Reference, options.RelativeOrientation2});
                self.relativeOrientation2Reference = options.RelativeOrientation2Reference;
            else
                self.relativeOrientation2 = options.RelativeOrientation2;
            end
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = generateMBDynInputString@mbdyn.pre.twoNodeJoint (self);
            
            str = self.addOutputLine (str, sprintf('%d', self.node1.label), 2, true, 'node 1 label');
            
            if ~isempty (self.relativeOffset1)
                out = self.makeCellIfNot (self.relativeOffset1);
                str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 3, true);
            end
            
            if ~isempty (self.relativeOrientation1)
                out = self.makeCellIfNot (self.relativeOrientation1);
                str = self.addOutputLine (str, self.commaSepList ('orientation', out{:}), 3, true);
            end
            
            addcomma = ~isempty (self.relativeOffset2) || ~isempty (self.relativeOrientation2);
            str = self.addOutputLine (str, sprintf('%d', self.node2.label), 2, addcomma, 'node 2 label');
            
            if ~isempty (self.relativeOffset2)
                addcomma = ~isempty (self.relativeOrientation2);
                out = self.makeCellIfNot (self.relativeOffset2);
                str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 3, addcomma);
            end
            
            if ~isempty (self.relativeOrientation2)
                out = self.makeCellIfNot (self.relativeOrientation2);
                str = self.addOutputLine (str, self.commaSepList ('orientation', out{:}), 3, false);
            end
            
            str = self.addOutputLine (str, ';', 1, false, sprintf ('end %s', self.type));
            
            str = self.addRegularization (str);
            
        end
        
        function draw (self, varargin)
            
%             options.AxesHandle = self.drawAxesH;
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
            
%             ref_node = mbdyn.pre.reference (self.node1.absolutePosition, ...
%                                             self.node1.absoluteOrientation, ...
%                                             [], []);
%                                         
%             ref_joint = mbdyn.pre.reference (self.relativeOffset1, self.relativePositionOrientation1, [], [], 'Parent', ref_node);
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
        end
        
    end
    
end