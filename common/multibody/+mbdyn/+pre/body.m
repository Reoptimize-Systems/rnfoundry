classdef body < mbdyn.pre.element
    
    properties (GetAccess = public, SetAccess = protected)
        
        mass;                 % the mass of the body
        relativeCentreOfMass; % the location of the center of mass with respect to the node in the reference frame of the node
        inertiaMatrix;        % Inertia_matrix referred to the center of mass of the mass
        inertialOrientation; 
        nMasses;
        
        nodeAttached;
        
    end
    
    methods
        
        function self = body (mass, cog, inertiamat, node, varargin)
        
            options.InertialOrientation = [];
            options.STLFile = '';
            options.UseSTLName = false;
            
            options = parse_pv_pairs (options, varargin);
            
            % call superclass constructor
            self = self@mbdyn.pre.element ('STLFile', options.STLFile, ...
                                           'UseSTLName', options.UseSTLName);
            
            if ~(isscalar (mass) && isnumeric (mass))
                error ('mass should be a numeric scalar value');
            end
            
            self.checkCOGVector (cog, true);
            self.checkInertiaMatrix (inertiamat, true);
            self.checkIsStructuralNode (node, true);
            
            self.mass = mass;
            self.nodeAttached = node;
            if isempty (cog)
                self.relativeCentreOfMass = 'null';
            else
                self.relativeCentreOfMass = cog;
            end
            self.inertiaMatrix = inertiamat;
            
            if ~isempty (options.InertialOrientation)
                if ischar (options.InertialOrientation) ...
                    if ~strcmp (options.InertialOrientation, 'node')
                        error ('InertialOrientation must be an orientation matrix or the keyword ''node''')
                    end
                else
                    self.checkOrientationMatrix (options.InertialOrientation, true);
                end
            end
                
            self.inertialOrientation = self.getOrientationMatrix (options.InertialOrientation);
            
        end
        
        function str = generateOutputString (self)
            
            str = self.addOutputLine ('' , '', 1, false, 'one-mass body');
            
            % delete newline character and space from start
            str(1:2) = [];
            
            str = self.addOutputLine (str, sprintf('body : %d, %d', self.label, self.nodeAttached.label), 1, true, 'label, node label');
            
            str = self.addOutputLine (str, self.commaSepList (self.mass), 2, true, 'mass');
            
            str = self.addOutputLine (str, self.commaSepList (self.relativeCentreOfMass), 2, true, 'relative centre of mass');
            
            addcomma = ~isempty (self.inertialOrientation);
            str = self.addOutputLine (str, self.commaSepList (self.inertiaMatrix), 2, addcomma, 'inertia matrix');
            
            if ~isempty (self.inertialOrientation)
                str = self.addOutputLine (str, self.commaSepList ('inertial', self.inertialOrientation), 2, false);
            end
            
            str = self.addOutputLine (str, ';', 1, false, 'end one-mass body');
            
        end
        
        function hax = draw (self, varargin)
            
            options.AxesHandle = [];
            options.ForceRedraw = false;
            options.Mode = 'solid';
            options.Light = false;
            
            options = parse_pv_pairs (options, varargin);
            
            hax = draw@mbdyn.pre.element ( self, ...
                    'AxesHandle', options.AxesHandle, ...
                    'ForceRedraw', options.ForceRedraw, ...
                    'Mode', options.Mode, ...
                    'Light', options.Light );

            self.setTransform ();
            
        end
        
    end
    
    methods (Access = protected)
        
        function setTransform (self)
            
            ref_node = mbdyn.pre.reference (self.nodeAttached.absolutePosition, ...
                                            self.nodeAttached.absoluteOrientation, ...
                                            [], []);
                                        
            ref_cog = mbdyn.pre.reference (self.relativeCentreOfMass, [], [], [], 'Parent', ref_node);
            
            M = [ ref_cog.orientm.orientationMatrix , ref_cog.pos; ...
                  0, 0, 0, 1 ];
            
            % matlab uses different convention to mbdyn for rotation
            % matrix
            M = self.mbdynOrient2Matlab (M);
                  
            set ( self.transformObject, 'Matrix', M );
            
        end        
        
    end
    
end