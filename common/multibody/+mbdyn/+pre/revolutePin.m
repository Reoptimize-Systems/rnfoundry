classdef revolutePin < mbdyn.pre.singleNodeJoint
    
    
    properties (GetAccess = public, SetAccess = protected)
        
        relativeOffset;
        nodeRelativeOrientation;
            
        pinPosition;
        absolutePinOrientation;
        
        initialTheta;
    end
    
    
    methods
        
        function self = revolutePin (node, relative_offset, absolute_pin_position, varargin)
            
            options.NodeRelativeOrientation = [];
            options.PinOrientation =  [];
            options.InitialTheta = [];
            
            options = parse_pv_pairs (options, varargin);
            
            % call the superclass constructor
            self = self@mbdyn.pre.singleNodeJoint (node);
            
            self.checkCartesianVector (relative_offset, true);
            self.checkCartesianVector (absolute_pin_position, true);
            self.checkOrientationMatrix (options.NodeRelativeOrientation, true);
            self.checkOrientationMatrix (options.PinOrientation, true);

            self.relativeOffset = relative_offset;
            self.nodeRelativeOrientation = self.getOrientationMatrix (options.NodeRelativeOrientation);

            self.pinPosition = absolute_pin_position;
            self.absolutePinOrientation = self.getOrientationMatrix (options.PinOrientation);
            
            self.initialTheta = options.InitialTheta;
            
        end
        
        function str = generateOutputString (self)
            
            str = self.addOutputLine ('' , '', 1, false, 'revolute pin');
            
            % delete newline character and space from start
            str(1:2) = [];
            
            str = self.addOutputLine (str , 'revolute pin', 1, true);
            
            str = self.addOutputLine (str, sprintf('%d', self.node.label), 2, true, 'node label');
            
            out = self.makeCellIfNot (self.relativeOffset);
            str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 3, true, 'node relative position' );
            
            if ~isempty (self.nodeRelativeOrientation)
                out = self.makeCellIfNot (self.nodeRelativeOrientation);
                str = self.addOutputLine (str, self.commaSepList ('orientation', out{:}), 3, true, 'node relative orientation');
            end
            
            out = self.makeCellIfNot (self.pinPosition);
            addcomma = ~(isempty (self.absolutePinOrientation) && isempty (self.initialTheta));
            str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 2, addcomma, 'pin absolute position');
            
            if ~isempty (self.absolutePinOrientation)
                addcomma = ~isempty (self.initialTheta);
                out = self.makeCellIfNot (self.absolutePinOrientation);
                str = self.addOutputLine (str, self.commaSepList ('orientation', out{:}), 2, addcomma, 'pin absolute orientation');
            end
            
            if ~isempty (self.initialTheta)
                out = self.makeCellIfNot (self.initialTheta);
                str = self.addOutputLine (str, self.commaSepList ('initial theta', out{:}), 2, false, 'initial theta');
            end
            
            str = self.addOutputLine (str, ';', 1, false, 'end revolute pin');
            
        end
        
    end
    
end