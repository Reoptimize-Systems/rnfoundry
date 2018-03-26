classdef genel < mbdyn.pre.element
    
    properties (GetAccess = public, SetAccess = protected)
        
        
    end
    
    methods
        
        function self = genel ()
        
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = sprintf ('    genel : %d, %s,', self.label, self.type);
            
        end
        
        function hax = draw (self, varargin)
            
            options.AxesHandle = [];
            options.ForceRedraw = false;
            options.Mode = 'solid';
            options.Light = false;
            
            options = parse_pv_pairs (options, varargin);
            
            hax = draw@mbdyn.pre.element ( self, ...
                    'AxesHandle', options.AxesHandle, ...
                    'ForceRedraw', options.ForceRedraw, ...
                    'Mode', options.Mode, ...
                    'Light', options.Light );

            self.setTransform ();
            
        end
        
    end
    
    methods (Access = protected)
        
        
        
    end
    
end