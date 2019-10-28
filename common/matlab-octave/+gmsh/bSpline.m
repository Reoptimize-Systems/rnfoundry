classdef bSpline < gmsh.curve
    
    properties
        
        controlPoints;
        
    end
    
    methods
        
        function self = bSpline (points, varargin)

            assert (isa (points, 'gmsh.point') && numel (points) >= 2, 'points must be an array of two or more gmsh.point objects');
            
            self.controlPoints = points;
            
            self.geoName = 'BSpline';
            
        end
        
        function str = generateGeoFileStr (self)
            
            str = generateGeoFileStr@gmsh.curve (self);
            
            point_tags = self.formatInteger (self.controlPoints.tag);
            
            str = sprintf ('%s (%s)', ...
                            str, ...
                            self.commaSepList (point_tags{:}) );
            
        end        
        
    end
    
end