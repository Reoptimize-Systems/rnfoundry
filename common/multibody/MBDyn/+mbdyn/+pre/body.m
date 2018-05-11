classdef body < mbdyn.pre.element
    
    properties (GetAccess = public, SetAccess = protected)
        
        mass;                 % the mass of the body
        relativeCentreOfMass; % the location of the center of mass with respect to the node in the reference frame of the node
        inertiaMatrix;        % Inertia_matrix referred to the center of mass of the mass
        inertialOrientation; 
        nMasses;
        
        nodeAttached;
        
    end
    
    properties (GetAccess = protected, SetAccess = protected)
        isPointMass;
    end
    
    methods
        
        function self = body (mass, cog, inertiamat, node, varargin)
            % constructs a single mass lumped rigid body or point mass
            %
            % Syntax
            %
            % bd = mbdyn.pre.body (mass, cog, inertiamat, node)
            % bd = mbdyn.pre.body (..., 'Parameter', value)
            %
            % Description
            %
            % mbdyn.pre.body constructs a lumped rigid body when connected
            % to a regular, 6 degree of freedom structural node, or a point
            % mass when connected to a rotationless, 3 degree of freedom
            % structural node.
            %
            % Input
            %
            %  mass - mass of the body
            %
            %  cog - (3 x 1) matrix containing the location of the centre
            %    of mass of the body relative to the attached node.
            %
            %  inertiamat - Inertia_matrix referred to the center of mass
            %    of the body. can be a (3 x 3) matrix or
            %    mbdyn.pre.orientmat object. The inertia_matrix is always
            %    referred to the center of mass of the body. However, it
            %    can be rotated locally using the 'InertialOrientation'
            %    options (see below.
            %
            %  node - node to which the body is attached, can be either an
            %    mbdyn.pre.structuralNode6dof object or a
            %    mbdyn.pre.structuralNode3dof object.
            %
            % Additional arguments may be supplied as parameter-value
            % pairs.
            %
            %  'InertialOrientation' - string, or (3 x 3) matrix or
            %    mbdyn.pre.orientmat onject. This options is used to rotate
            %    the inertial matrix supplied in inertiamat. If this is a
            %    string, it must be 'node' which is the default, and
            %    assumes the inertia matrix is input in the node reference
            %    frame. Otherwise it is rotated according to the supplied
            %    inertial matrix.
            %
            %  'STLFile' - path to an STL file used to plot the body in
            %    visualisations. If not supplied, a default shape is used.
            %
            %  'UseSTLName' - true/false flag. If an stl file is loaded,
            %    this option determines whether the name stored in the
            %    stl file is assigned to the body. Default is false.
            %
            % Output
            %
            %  bd - mbdyn.pre.body object
            %
            %
            % See Also: mbdyn.pre.bodyMultiMass
            %

            [options, nopass_list] = mbdyn.pre.body.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            % call superclass constructor
            self = self@mbdyn.pre.element ( pvpairs{:} );
            
            if ~(isscalar (mass) && isnumeric (mass))
                error ('mass should be a numeric scalar value');
            end
            
            self.isPointMass = false;
            if isa (node, 'mbdyn.pre.structuralNode3dof')
               self.isPointMass = true; 
            end
            
            self.mass = mass;
            self.nodeAttached = node;
            
            if ~self.isPointMass
                
                self.checkCOGVector (cog, true);
                self.checkInertiaMatrix (inertiamat, true);
                self.checkIsStructuralNode (node, true);
            
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
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for a body
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (bd)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  bd - mbdyn.pre.body object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = self.addOutputLine ('' , '', 1, false, 'one-mass body');
            
            % delete newline character and space from start
            str(1:2) = [];
            
            str = self.addOutputLine (str, sprintf('body : %d, %d', self.label, self.nodeAttached.label), 1, true, 'label, node label');

            addcomma = ~self.isPointMass;
            str = self.addOutputLine (str, self.commaSepList (self.mass), 2, addcomma, 'mass');
                
            if ~self.isPointMass

                str = self.addOutputLine (str, self.commaSepList (self.relativeCentreOfMass), 2, true, 'relative centre of mass');

                addcomma = ~isempty (self.inertialOrientation);
                str = self.addOutputLine (str, self.commaSepList (self.inertiaMatrix), 2, addcomma, 'inertia matrix');

                if ~isempty (self.inertialOrientation)
                    str = self.addOutputLine (str, self.commaSepList ('inertial', self.inertialOrientation), 2, false);
                end
            
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
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.element.defaultConstructorOptions ();
            
            parentfnames = fieldnames (options);
            
            options.InertialOrientation = [];
            
            allfnames = fieldnames (options);
            
            nopass_list = setdiff (allfnames, parentfnames);
            
        end
        
    end
    
end