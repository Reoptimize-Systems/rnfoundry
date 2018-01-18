classdef clamp < mbdyn.pre.singleNodeJoint
    
    
    properties (GetAccess = public, SetAccess = protected)
        
        position;
        orientation;
            
        positionReference;
        orientationReference;
        
    end
    
    properties (Dependent)
        
        absolutePosition;
        absoluteOrientation;
        
    end
    
    properties (GetAccess = protected, SetAccess = protected)
        posIsNode;
        orientIsNode;
    end
    
    methods
        
        function self = clamp (node, position, orientation, varargin)
            
            options.PositionReference = 'node';
            options.OrientationReference = 'node';
            
            options = parse_pv_pairs (options, varargin);
            
            % call the superclass constructor
            self = self@mbdyn.pre.singleNodeJoint (node);
            
            self.type = 'clamp';
            self.posIsNode = false;
            self.orientIsNode = false;
            
            if ~isempty (position)
                if ischar (position)
                    if strcmp (position, 'node')
                        self.position = position;
                        self.posIsNode = true;
                    else
                        error ('Clamp position must be the keyword ''node'' or a position vector');
                    end
                else
                    self.checkCartesianVector (position, true);
                    self.checkNodeReferenceType (options.PositionReference, true);
                    self.position = {'reference', options.PositionReference, options.Position};
                end
            else
                self.position = [];
            end
            
            if ~isempty (orientation)
                if ischar (orientation)
                    if strcmp (orientation, 'node')
                        self.orientation = orientation;
                        self.orientIsNode = true;
                    else
                        error ('Clamp orientation must be the keyword ''node'' or an orientation matrix or mbdyn.pre.orientm object');
                    end
                else
                    self.checkOrientationMatrix (orientation, true);
                    self.checkNodeReferenceType (options.OrientationReference, true);
                    self.orientation = {'reference', options.OrientationReference, self.getOrientationMatrix(options.Orientation)};
                end
            else
                self.orientation = [];
            end
            
        end
        
        function str = generateOutputString (self)
            
            str = generateOutputString@mbdyn.pre.singleNodeJoint(self);
            
            addcomma = ~isempty (self.position);
            str = self.addOutputLine (str, sprintf('%d', self.node.label), 2, addcomma, 'node label');
            
            if ~isempty (self.position)
                addcomma = ~isempty (self.orientation);
                out = self.makeCellIfNot (self.position);
                str = self.addOutputLine (str, self.commaSepList ('position', out{:}), 3, addcomma, 'absolute position');
            end
            
            if ~isempty (self.orientation)
                out = self.makeCellIfNot (self.orientation);
                str = self.addOutputLine (str, self.commaSepList ('orientation', out{:}), 3, false, 'absolute orientation');
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
        
        function abspos = get.absolutePosition (self)
            % gets the position of the clamp in the global frame
            
            if self.posIsNode
                    
                abspos = self.node.absolutePosition;
                
            else
                
                switch self.positionReference

                    case 'global'

                        abspos = self.position;

                    case {'node', 'local'}

                        abspos = self.node.relativeToAbsolutePosition (self.position);

                end
            
            end
            
        end
        
        function absorientm = get.absoluteOrientation (self)
            % gets the orientation of the clamp in the global frame
            
            if self.orientIsNode
                
                absorientm = self.node.absoluteOrientation;
                
            else
                switch self.orientationReference

                    case 'global'

                        absorientm = self.orientation;

                    case {'node', 'local'}

                        absorientm = self.node.relativeToAbsoluteOrientation (self.orientation);

                end
            end
            
        end
        
    end
    
    methods (Access = protected)
        
        function setTransform (self)
            
            M = [ self.absoluteOrientation.orientationMatrix, self.absolutePosition; ...
                  0, 0, 0, 1 ];
                  
            set ( self.transformObject, 'Matrix', M );
            
        end
        
    end
    
    
    
end