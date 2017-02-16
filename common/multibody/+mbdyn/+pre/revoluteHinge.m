classdef revoluteHinge < mbdyn.pre.twoNodeJoint
    
    
    properties (GetAccess = public, SetAccess = protected)
        
        relativeOffset1;
        relativeOrientation1;
            
        relativeOffset2;
        relativeOrientation2;
        
    end
    
    
    methods
        
        function self = revoluteHinge (node1, node2, position1, position2, varargin)
            
            options.RelativeOrientation1 =  [];
            options.RelativeOrientation2 =  [];
            options.InitialTheta = [];
            options.Friction = [];
            options.Preload = [];
            options.FrictionModel = [];
            options.ShapeFunction = [];
            
            options = parse_pv_pairs (options, varargin);
            
            % call the superclass constructor
            self = self@mbdyn.pre.twoNodeJoint (node1, node2);
            
            self.checkJointPositionOffset (position1);
            self.checkJointPositionOffset (position2);
            self.checkJointOrientation (options.RelativeOrientation1);
            self.checkJointOrientation (options.RelativeOrientation2);
            
            self.relativeOffset1 = position1;
            self.relativeOrientation1 = self.getOrientationMatrix (options.RelativeOrientation1);

            self.relativeOffset2 = position2;
            self.relativeOrientation2 = self.getOrientationMatrix (options.RelativeOrientation2);
            
        end
        
        function str = generateOutputString (self)
            
            str = self.addOutputLine ('' , '', 1, false, 'revolute hinge');
            
            % delete newline character and space from start
            str(1:2) = [];
            
            str = self.addOutputLine (str , 'revolute hinge', 1, true);
            
            str = self.addOutputLine (str, sprintf('%d', self.node1.label), 2, true, 'node 1 label');
            
            out = self.makeCellIfNot (self.relativeOffset1);
            str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 3, true);
            
            if ~isempty (self.relativeOrientation1)
                out = self.makeCellIfNot (self.relativeOrientation1);
                str = self.addOutputLine (str, self.commaSepList ('orientation', out{:}), 3, true);
            end
            
            str = self.addOutputLine (str, sprintf('%d', self.node2.label), 2, true, 'node 2 label');
            
            out = self.makeCellIfNot (self.relativeOffset2);
            addcomma = ~isempty (self.relativeOrientation2);
            str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 3, addcomma);
            
            if ~isempty (self.relativeOrientation2)
                out = self.makeCellIfNot (self.relativeOrientation2);
                str = self.addOutputLine (str, self.commaSepList ('orientation', out{:}), 3, false);
            end
            
            str = self.addOutputLine (str, ';', 1, false, 'end total joint');
            
        end
        
    end
    
end