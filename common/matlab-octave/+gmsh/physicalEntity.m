classdef physicalEntity < gmsh.base
    
    properties
        
        label;
        
    end

    
    methods
        
        function self = physicalEntity (varargin)
            
            options.Label = '';
            
            options = parse_pv_pairs (options, varargin);
            
            if ~isempty (options.Label), check.isCharOrScalarString (options.Label, true, 'Label'); end
            
            self.label = options.Label;
            
        end
        
        function str = generateGeoFileStr (self)
            
            if isempty (self.label)
                str = sprintf ('%s (%s) =', self.geoName, self.formatInteger (self.tag));
            else
                str = sprintf ('%s ("%s", %s) =', self.geoName, self.label, self.formatInteger (self.tag));
            end
            
        end
        
    end
    
end