classdef line < gmsh.curve
    
    properties
        
        point1;
        point2;
        
    end
    
    methods
        
        function self = line (p1, p2, varargin)

            assert (isa (p1, 'gmsh.point') && isscalar (p1), 'p1 must be a scalar gmsh.point object');
            assert (isa (p2, 'gmsh.point') && isscalar (p2), 'p2 must be a scalar gmsh.point object');
            
            self.point1 = p1;
            self.point2 = p2;
            
            self.geoName = 'Line';
            
        end
        
        function str = generateGeoFileStr (self)
            
            str = generateGeoFileStr@gmsh.curve (self);
            
            str = sprintf ('%s (%s, %s)', ...
                            str, ...
                            self.formatInteger (self.point1.tag), ...
                            self.formatInteger (self.point2.tag) );
            
        end        
        
    end
    
end