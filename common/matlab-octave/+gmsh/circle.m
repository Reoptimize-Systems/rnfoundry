classdef circle < gmsh.curve
    
    properties
        
        startPoint;
        centrePoint;
        endPoint;
        
    end
    
    methods
        
        function self = circle (start_point, centre_point, end_point, varargin)

            assert (isa (start_point, 'gmsh.point') && isscalar (start_point), 'start_point must be a scalar gmsh.point object');
            assert (isa (centre_point, 'gmsh.point') && isscalar (centre_point), 'centre_point must be a scalar gmsh.point object');
            assert (isa (end_point, 'gmsh.point') && isscalar (end_point), 'end_point must be a scalar gmsh.point object');
            
            self.startPoint = start_point;
            self.centrePoint = centre_point;
            self.endPoint = end_point;
            
            self.geoName = 'Circle';
            
        end
        
        function str = generateGeoFileStr (self)
            
            str = generateGeoFileStr@gmsh.curve (self);

            str = sprintf ('%s (%s, %s, %s)', ...
                            str, ...
                            self.formatInteger (self.startPoint.tag), ...
                            self.formatInteger (self.centrePoint.tag), ...
                            self.formatInteger (self.endPoint.tag) );
            
        end
        
    end
    
end