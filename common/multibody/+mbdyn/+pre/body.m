classdef body < mbdyn.pre.element
    
    properties (GetAccess = public, SetAccess = protected)
        
        mass;
        relativeCentreOfMass;
        inertiaMatrix;
        inertialOrientation;
        
        nodeAttached;
        
    end
    
    methods
        
        function self = body (mass, cog, inertiamat, node, varargin)
        
            options.InertialOrientation = [];
            
            options = parse_pv_pairs (options, varargin);
            
            if ~(isscalar (mass) && isnumeric (mass))
                error ('mass should be a numeric scalar value');
            end
            
            self.checkCOGVector (cog, true);
            self.checkInertiaMatrix (inertiamat, true);
            self.checkIsStructuralNode (node, true);
            
            self.mass = mass;
            if isempty (cog)
                self.relativeCentreOfMass = 'null';
            else
                self.relativeCentreOfMass = cog;
            end
            self.inertiaMatrix = inertiamat;
            self.nodeAttached = node;
            
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
        
    end
    
    methods (Access = protected)
        
        function ok = checkInertiaMatrix (self, mat, throw)
            
            ok = self.check3X3Matrix (mat, false);
            
            if ~ok && throw
                error ('Inertia matrix must be a 3 x 3 numeric matrix');
            end

        end
        
        function ok = checkCOGVector (self, cog, throw)

            if isempty (cog) || (ischar (cog) && strcmp (cog, 'null'))
                ok = true;
            else
                ok = self.checkCartesianVector (cog, false);
            end
            
            if ~ok && throw
                error ('Centre of gravity offset must 3 element numeric column vector, or keyword ''null'' or empty');
            end
                
            
        end
        
        
    end
    
end