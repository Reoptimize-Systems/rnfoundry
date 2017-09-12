classdef bodyMultiMass < mbdyn.pre.element
    
    properties (GetAccess = public, SetAccess = protected)
        
        mass;                 % the masses of the bodies
        relativeCentreOfMass; % the location of the centers of mass with respect to the node in the reference frame of the node
        inertiaMatrix;        % Inertia_matrices referred to the center of mass of each mass
        inertialOrientation;  
        nMasses;
        
        nodeAttached;
        
    end
    
    properties (GetAccess = protected, SetAccess = protected)
        isPointMass;
        bodies;
    end
    
    methods
        
        function self = bodyMultiMass (mass, cog, inertiamat, node, varargin)
        
            options.InertialOrientation = {};
            options.STLFiles = {};
            options.UseSTLName = false;
            
            options = parse_pv_pairs (options, varargin);
            
            % call superclass constructor
%             self = self@mbdyn.pre.element ('STLFile', options.STLFile, ...
%                                            'UseSTLName', options.UseSTLName);
            
            if ~( isnumeric (mass) && isreal (mass) )
                error ('mass should be a numeric array of one or more mass values.');
            end
            

            if ~ ( iscell (cog) ...
                    && iscell (inertiamat) ...
                    && samesize (mass, cog, inertiamat) )
                error ('cog must be a cell array of 3 element vectors of the same size as mass, and inertiamat must a cell array of inertia matrices.');
            end

            self.checkIsStructuralNode (node, true);
            
            self.isPointMass = false;
            if isa (node, 'mbdyn.pre.structuralNode3dof')
               self.isPointMass = true; 
            end
            
            for ind = 1:numel (mass)
                
                if ~self.isPointMass
                    self.checkCOGVector (cog{ind}, true);
                    self.checkInertiaMatrix (inertiamat{ind}, true);

                    if isempty (cog{ind})
                        self.relativeCentreOfMass{ind} = 'null';
                    else
                        self.relativeCentreOfMass{ind} = cog{ind};
                    end
                
                    self.inertiaMatrix{ind} = inertiamat{ind};
                    
                    
                end
            
            end
            
            self.mass = mass;
            self.nodeAttached = node;
            
            if isempty (options.InertialOrientation)
                options.InertialOrientation = cell (size(mass));
            else
                if ~( iscell (options.InertialOrientation) ...
                        && samesize (mass, options.InertialOrientation) )
                    error ('InertialOrientation has been supplied, but is not a cell array the same size as mass');
                end
            end
                
            for ind = 1:numel (self.mass)
                if ~isempty (options.InertialOrientation{ind})
                    if ischar (options.InertialOrientation{ind}) ...
                        if ~strcmp (options.InertialOrientation{ind}, 'node')
                            error ('InertialOrientation must be an orientation matrix or the keyword ''node''')
                        end
                    else
                        self.checkOrientationMatrix (options.InertialOrientation{ind}, true);
                    end
                end
                self.inertialOrientation{ind} = self.getOrientationMatrix (options.InertialOrientation{ind});
            end

            if isempty (options.STLFiles)
                options.STLFiles = cell (size(mass));
            else
                if ~ (iscell (options.STLFiles) && samesize (options.STLFiles, self.masses) )
                    error ('STLFiles must be a cell array of the smae size as the mass.');
                end
            end
            
            % create array of bodies for drawing
            self.bodies = mbdyn.pre.body ( mass(1), cog{1}, inertiamat{1}, node, ...
                                        'STLFile', options.STLFiles{1}, ...
                                        'UseSTLName', options.UseSTLName, ...
                                        'InertialOrientation', self.inertialOrientation{1} );
            for ind = 2:numel (self.mass)
                self.bodies(ind) = mbdyn.pre.body ( mass(ind), cog{ind}, inertiamat{ind}, node, ...
                                        'STLFile', options.STLFiles{ind}, ...
                                        'UseSTLName', options.UseSTLName, ...
                                        'InertialOrientation', self.inertialOrientation{ind} );
            end
            
        end
        
        function str = generateOutputString (self)
            
                
            str = self.addOutputLine ('' , '', 1, false, 'multiple-mass body');

            % delete newline character and space from start
            str(1:2) = [];

            str = self.addOutputLine (str, sprintf('body : %d, %d', self.label, self.nodeAttached.label), 1, true, 'label, node label');

            str = self.addOutputLine (str, self.commaSepList ('condense'), 2, true);

            for ind = 1:numel (self.mass)

                str = self.addOutputLine (str, self.commaSepList (self.mass(ind)), 2, true, 'mass');

                str = self.addOutputLine (str, self.commaSepList (self.relativeCentreOfMass{ind}), 2, true, 'relative centre of mass');

                addcomma = ( ~isempty (self.inertialOrientation{ind}) || ind < numel (self.mass));
                str = self.addOutputLine (str, self.commaSepList (self.inertiaMatrix{ind}), 2, addcomma, 'inertia matrix');

                if ~isempty (self.inertialOrientation{ind})
                    str = self.addOutputLine (str, self.commaSepList ('inertial', self.inertialOrientation{ind}), 2, false);
                end

            end

            str = self.addOutputLine (str, ';', 1, false, 'end multiple-mass body');
            
        end
        
        function setSize (self, n, sx, sy, sz)
            % set size of the n'th mass for drawing
            
            self.bodies(n).setSize(sx, sy, sz);
        end
        
        function hax = draw (self, varargin)
            
            options.AxesHandle = [];
            options.ForceRedraw = false;
            options.Mode = 'solid';
            options.Light = false;
            
            options = parse_pv_pairs (options, varargin);
            
            hax = options.AxesHandle;
            
            for ind = 1:numel (self.bodies)
                
                hax = self.bodies(ind).draw ( ...
                        'AxesHandle', hax, ...
                        'ForceRedraw', options.ForceRedraw, ...
                        'Mode', options.Mode, ...
                        'Light', options.Light );
                    
                hold (hax, 'on');

            end
            
            hold (hax, 'off');
            
        end
        
    end
    
    
end