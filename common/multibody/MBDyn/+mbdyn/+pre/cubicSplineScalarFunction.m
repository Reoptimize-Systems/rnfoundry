classdef cubicSplineScalarFunction < mbdyn.pre.scalarFunction
    
    properties (GetAccess = public, SetAccess = protected)
        x;
        y;
        extrapolate;
    end
    
    methods
        
        function self = cubicSplineScalarFunction (name, x, y, varargin)
            % mbdyn.pre.cubicSplineScalarFunction constructor
            %
            %
            % Syntax
            %
            % cssf = mbdyn.pre.cubicSplineScalarFunction (name, c)
            %
            % Description
            %
            % Implements a scalar function cubic natural spline
            % interpolation between a set of x and y points. There must be
            % at least 3 points.
            %
            % Input
            %
            %  name - unique name for the function
            %
            %  x - scalar value of the function
            %
            %  y - scalar value of the function
            %
            % Output
            %
            %  cssf - mbdyn.pre.cubicSplineScalarFunction object
            %
            %
            % See Also: mbdyn.pre.scalarFunction
            %
            
            [ options, nopass_list ] = mbdyn.pre.cubicSplineScalarFunction.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs ( options, nopass_list);
            
            self = self@mbdyn.pre.scalarFunction (name, 'cubicspline');
            
            assert (isnumeric (x) && isreal (x) && isvector (x), ...
                    'x must be a vector of real numeric values' );
            assert (isnumeric (x) && isreal (x) && isvector (x), ...
                    'y must be a vector of real numeric values' );
            assert (numel (x) == numel (y), 'x must be the same length as y');
            assert (numel(x) >= 3, 'There must be at least 3 points for the interpolation');
            
            self.checkLogicalScalar (options.Extrapolate, true, 'Extrapolate');
            
            self.name = name;
            self.extrapolate = options.Extrapolate;
            
            % store x and y ensuring they are row vectors as this makes it
            % easier to merge them for output as a comma separated list
            % later
            self.x = x(:)';
            self.y = y(:)';
            
        end
        
        function str = generateMBDynInputString (self)
            
            points = [self.x; self.y];
            
            points = points(:)';
            
            things = { ['"', self.name, '"'] };
            
            if self.extrapolate == false
                things = [ things, { 'do not extrapolate' } ];
            end
            
            things = [ things, { self.fcnType, ...
                                 points, ...
                                 'end' } ];
                             
            str = self.commaSepList ( things{:} );
            
        end
        
    end
    
    methods (Static)
        
        function [ options, nopass_list ] = defaultConstructorOptions ()
            
            options = mbdyn.pre.scalarFunction.defaultConstructorOptions ();
            
            parentfnames = fieldnames (options);
            
            % add default options common to all revoluteHinge objects
            options.Extrapolate = true;
            
            allfnames = fieldnames (options);
            
            C = setdiff (allfnames, parentfnames, 'stable');
            
            nopass_list = C;
            
        end
        
    end
    
end