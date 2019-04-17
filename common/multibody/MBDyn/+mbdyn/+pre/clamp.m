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
            % clamp joint constructor
            %
            % Syntax
            %
            % cl = clamp (node, position, orientation)
            % cl = clamp (node, position, orientation, 'Parameter', value)
            %
            % Description
            %
            % mbdyn.pre.clamp is a joint which grounds all 6 degrees of
            % freedom of a node in an arbitrary position and orientation
            % that remains fixed.
            %
            % Input
            %
            %  node - mbdyn.pre.structuralNode object representing the node
            %   to be clamped
            %
            %  position - either a string 'node', or a (3 x 1) vector. The
            %   keyword 'node' forces the joint to use the node’s position.
            %   Otherwise the vector contains the position in the absolute
            %   frame in which the node is to be clamped. To specify an
            %   alternative reference frame the PositionReference option
            %   can be used (see below).
            %
            %  orientation - either a string 'node', or an
            %   mbdyn.pre.orientmat object. The keyword 'node' forces the
            %   joint to use the node’s orientation. Otherwise the
            %   mbdyn.pre.orientmat object represents the orientation in
            %   the absolute frame in which the node is to be clamped. To
            %   specify an alternative reference frame the
            %   OrientationReference option can be used (see below).
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'PositionReference' - optional string containing an
            %    alternative reference for te clamp position. Valid strings
            %    are: 'node', 'local' and 'global'.
            %
            %  'OrientationReference' - optional string containing an
            %    alternative reference for the clamp orientation. Valid
            %    strings are: 'node', 'local' and 'global'.
            %
            % Output
            %
            %  cl - mbdyn.pre.clamp object
            %
            %
            %
            % See Also: 
            %
            
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
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for clamp joint
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (cl)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  cl - mbdyn.pre.clamp object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = generateMBDynInputString@mbdyn.pre.singleNodeJoint(self);
            
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