classdef point < gmsh.base
    
    properties
        
        x;
        y;
        z;
        elSize;
        
    end
    
    methods
        
        function self = point (x, y, z, varargin)
            
            options.ElementSize = [];
            
            options = parse_pv_pairs (options, varargin);
            
            check.isNumericScalar (x, true, 'x');
            check.isNumericScalar (y, true, 'y');
            check.isNumericScalar (z, true, 'z');
            if ~isempty (options.ElementSize),  check.isNumericScalar (x, true, 'x'); end
            
            self.x = x;
            self.y = y;
            self.z = z;
            self.elSize = options.ElementSize;
            self.geoName = 'Point';
            
        end
        
        function str = generateGeoFileStr (self)
            
            str = generateGeoFileStr@gmsh.base (self);
            
            str = sprintf ('%s (%s, %s, %s', ...
                            str, ...
                            self.formatNumber (self.x), ...
                            self.formatNumber (self.y), ...
                            self.formatNumber (self.z) );
            
            if isempty (self.elSize)
                str = sprintf ('%s)', str);
            else
                str = sprintf ('%s, %s)', ...
                                str, ...
                                self.formatNumber (self.elSize) );
            end
            
        end        
        
    end
    
    
end