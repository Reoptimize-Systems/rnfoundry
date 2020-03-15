classdef linearScalarFunction < mbdyn.pre.scalarFunction
    
    properties (GetAccess = public, SetAccess = protected)
        x1;
        y1;
        x2;
        y2;
    end
    
    methods
        
        function self = linearScalarFunction (name, x1, y1, x2, y2, varargin)
            % mbdyn.pre.linearScalarFunction constructor
            %
            %
            % Syntax
            %
            % lsf = mbdyn.pre.linearScalarFunction (name, x1, y1, x2, y2)
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
            %  x1 - x coordinate of the first point
            %
            %  y1 - y coordinate of the first point
            %
            %  x2 - x coordinate of the second point
            %
            %  y2 - y coordinate of the second point
            %
            % Output
            %
            %  lsf - mbdyn.pre.linearScalarFunction object
            %
            %
            % See Also: mbdyn.pre.scalarFunction
            %
            
            [ options, nopass_list ] = mbdyn.pre.linearScalarFunction.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs ( options, nopass_list);
            
            self = self@mbdyn.pre.scalarFunction (name, 'linear', pvpairs{:});
            
            self.checkNumericScalar (x1, true, 'x1');
            self.checkNumericScalar (y1, true, 'y1');
            self.checkNumericScalar (x2, true, 'x2');
            self.checkNumericScalar (y2, true, 'y2');
            
            % store x and y ensuring they are row vectors as this makes it
            % easier to merge them for output as a comma separated list
            % later
            self.x1 = x1;
            self.y1 = y1;
            self.x2 = x2;
            self.y2 = y2;
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = self.commaSepList ( ['"', self.name, '"'], ...
                                      self.fcnType, ...
                                      self.x1, ...
                                      self.y1, ...
                                      self.x2, ...
                                      self.y2 );
            
        end
        
    end
    
    methods (Static)
        
        function [ options, nopass_list ] = defaultConstructorOptions ()
            
            options = mbdyn.pre.scalarFunction.defaultConstructorOptions ();
            
            nopass_list = {};
            
        end
        
    end
    
end