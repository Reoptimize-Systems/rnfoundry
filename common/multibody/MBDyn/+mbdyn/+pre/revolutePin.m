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
            options.NodeOffsetReference = 'global';
            options.NodeOrientationReference = 'global';
            
            options = parse_pv_pairs (options, varargin);
            
            % call the superclass constructor
            self = self@mbdyn.pre.singleNodeJoint (node);
            
            self.type = 'revolute pin';
            
            self.checkCartesianVector (relative_offset, true);
            self.checkCartesianVector (absolute_pin_position, true);
            self.checkOrientationMatrix (options.NodeRelativeOrientation, true);
            self.checkOrientationMatrix (options.PinOrientation, true);
            self.checkNodeReferenceType (options.NodeOffsetReference, true);
            self.checkNodeReferenceType (options.NodeOrientationReference, true);

            self.relativeOffset = {'reference', options.NodeOffsetReference, relative_offset};
            self.nodeRelativeOrientation = {'reference', options.NodeOrientationReference, self.getOrientationMatrix(options.NodeRelativeOrientation)};

            self.pinPosition = absolute_pin_position;
            self.absolutePinOrientation = self.getOrientationMatrix (options.PinOrientation);
            
            self.initialTheta = options.InitialTheta;
            
        end
        
        function str = generateOutputString (self)
            
            str = generateOutputString@mbdyn.pre.singleNodeJoint(self);
            
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
            
            str = self.addOutputLine (str, ';', 1, false, sprintf('end %s', self.type));
            
        end
        
        function draw (self, varargin)
            
            options.AxesHandle = self.drawAxesH;
            options.ForceRedraw = false;
            options.Mode = 'solid';
            
            options = parse_pv_pairs (options, varargin);
            
            draw@mbdyn.pre.element ( self, ...
                'AxesHandle', options.AxesHandle, ...
                'ForceRedraw', options.ForceRedraw, ...
                'Mode', options.Mode );

            self.setTransform ();
            
        end
        
    end
    
    methods (Access = protected)
        
        function setTransform (self)
            
            ref_pin = mbdyn.pre.reference ( self.pinPosition, ...
                                            mbdyn.pre.orientmat ('orientation', self.absolutePinOrientation), ...
                                            [], []);
                                        
%             ref_joint = mbdyn.pre.reference (self.relativeOffset, ...
%                 mbdyn.pre.orientmat ('orientation', self.nodeRelativeOrientation), ...
%                 [], [], 'Parent', ref_pin);
            
            M = [ ref_pin.orientm.orientationMatrix , ref_pin.pos; ...
                  0, 0, 0, 1 ];
            
            % matlab uses different convention to mbdyn for rotation
            % matrix
            M = self.mbdynOrient2Matlab (M);
                  
            set ( self.transformObject, 'Matrix', M );
            
        end
        
    end
    
    
    
end