classdef joint < mbdyn.pre.element
    
    properties
        
    end
    
    methods
        
        function self = joint (varargin)
            % generic base class for joints which constrains two nodes
            
            [options, nopass_list] = mbdyn.pre.joint.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            % call superclass constructor
            self = self@mbdyn.pre.element ( pvpairs{:} );
            
        end
        
        function str = generateMBDynInputString (self)
            str = sprintf ('    joint : %d, %s,', self.label, self.type);
        end
        
    end
    
    methods (Access = protected)
        

        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.element.defaultConstructorOptions ();
            
            nopass_list = {};
            
        end
        
    end
    
end