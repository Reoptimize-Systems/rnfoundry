classdef spline < gmsh.curve
    
    properties
        
        points;
        
    end
    
    methods
        
        function self = spline (points, varargin)

            assert (isa (points, 'gmsh.point') && numel (points) >= 2, 'points must be an array of two or more gmsh.point objects');
            
            self.points = points;
            
            self.geoName = 'Spline';
            
        end
        
        function str = generateGeoFileStr (self)
            
            str = generateGeoFileStr@gmsh.curve (self);
            
            point_tags = self.formatInteger (self.points.tag);
            
            str = sprintf ('%s (%s)', ...
                            str, ...
                            self.commaSepList (point_tags{:}) );
            
        end        
        
    end
    
end