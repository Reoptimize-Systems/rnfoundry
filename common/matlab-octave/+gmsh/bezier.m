classdef bezier < gmsh.curve
    
    properties
        
        controlPoints;
        
    end
    
    methods
        
        function self = bezier (points, varargin)

            assert (isa (points, 'gmsh.point') && numel (points) >= 3, 'points must be an array of three or more gmsh.point objects');
            
            self.controlPoints = points;
            
            self.geoName = 'Bezier';
            
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