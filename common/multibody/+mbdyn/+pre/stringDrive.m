classdef stringDrive < mbdyn.pre.drive
    
    properties (GetAccess = public, SetAccess = protected)
        
        string;
        
    end
    
    methods
        
        function self = stringDrive (string, varargin)
            
            error ('not yet implemented')
        
            if ~ischar (string)
                error ('''string'' must be a char array')
            end
            
            self.string = string;
                
        end
        
    end
    
end