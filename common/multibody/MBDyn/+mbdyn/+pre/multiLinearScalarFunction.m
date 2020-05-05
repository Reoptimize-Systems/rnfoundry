classdef multiLinearScalarFunction < mbdyn.pre.scalarFunction
    
    properties (GetAccess = public, SetAccess = protected)
        x;
        y;
        extrapolate;
    end
    
    methods
        
        function self = multiLinearScalarFunction (name, x, y, varargin)
            % mbdyn.pre.multiLinearScalarFunction constructor
            %
            %
            % Syntax
            %
            % mlsf = mbdyn.pre.multiLinearScalarFunction (name, x, y)
            %
            % Description
            %
            % Implements a linear function defined by a line passing
            % through two points.
            %
            % Input
            %
            %  name - unique name for the function
            %
            %  x - vector of x coordinates of the points
            %
            %  y - vector of y coordinates of the points
            %
            % Output
            %
            %  mlsf - mbdyn.pre.multiLinearScalarFunction object
            %
            %
            % See Also: mbdyn.pre.scalarFunction
            %
            
            [ options, nopass_list ] = mbdyn.pre.multiLinearScalarFunction.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs ( options, nopass_list);
            
            self = self@mbdyn.pre.scalarFunction (name, 'multilinear', pvpairs{:});
            
            assert (isnumeric (x) && isreal(x), 'x must be a real numeric vector');
            assert (isnumeric (y) && isreal(y), 'y must be a real numeric vector');
            assert (numel (x) == numel (y), 'x must be the same length as y');
            mbdyn.pre.base.checkLogicalScalar (options.Extrapolate, true, 'Extrapolate');
            
            % store x and y ensuring they are row vectors as this makes it
            % easier to merge them for output as a comma separated list
            % later
            self.x = x(:).';
            self.y = y(:).';
            self.extrapolate = options.Extrapolate;
            
        end
        
        function str = generateMBDynInputString (self)
            
            tmp = [ self.x; self.y ];
            
            if self.extrapolate
            
                str = self.commaSepList ( ['"', self.name, '"'], ...
                                          self.fcnType, ...
                                          tmp(:), ... 
                                          'end' );
                                  
            else
                
                str = self.commaSepList ( ['"', self.name, '"'], ...
                                          self.fcnType, ...
                                          'do not extrapolate', ...
                                          tmp(:), ... 
                                          'end' );
                
            end
            
        end
        
    end
    
    methods (Static)
        
        function [ options, nopass_list ] = defaultConstructorOptions ()
            
            options = mbdyn.pre.scalarFunction.defaultConstructorOptions ();
            
            options.Extrapolate = true;
            
            nopass_list = {'Extrapolate'};
            
        end
        
    end
    
end