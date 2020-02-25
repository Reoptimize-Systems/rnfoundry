classdef drive < mbdyn.pre.base
    % generic base class for all mbdyn drives
    
    properties
        name;
    end
    
    methods
        
        function self = drive (varargin)
            % construct a generic mbdyn drive object
            
            [options, ~] = mbdyn.pre.drive.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            assert (ischar (options.Name), 'Name must be a character vector');
            
            self = self@mbdyn.pre.base ();
            
            self.name = options.Name;
            self.type = 'drive';
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = generateMBDynInputString@mbdyn.pre.element (self);
            
            str = sprintf ('%s    drive : %d, %s,', str, self.label, self.type);
            
        end
        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options.Name = '';
            
            nopass_list = {};
            
        end
        
    end
    
end