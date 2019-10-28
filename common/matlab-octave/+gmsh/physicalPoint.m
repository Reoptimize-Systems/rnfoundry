classdef physicalPoint < gmsh.physicalEntity
    
    properties
        
        points;
        
    end
    
    methods
        
        function self = physicalPoint (points, varargin)
            
            options.Label = '';
            
            options = parse_pv_pairs (options, varargin);
            
            
            self = self@gmsh.physicalEntity ('Label', options.Label);
            
            assert (isa (points, 'gmsh.point'), 'points must be an array of one or more gmsh.point objects');
            
            self.points = points;
            self.geoName = 'Physical Point';
            
        end
        
        function str = generateGeoFileStr (self)
            
            str = generateGeoFileStr@gmsh.physicalEntity (self);
            
            point_tags = self.formatInteger (self.points.tag);
            
            str = sprintf ('%s (%s)', ...
                            str, ...
                            self.commaSepList (point_tags{:}) );
            
        end
        
        
    end
    
    
end