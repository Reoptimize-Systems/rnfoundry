classdef tplDriveCaller < mbdyn.pre.driveCaller
    
    properties (GetAccess = public, SetAccess = protected)

        driveCallers;
        
    end
    
    methods
        
        function self = tplDriveCaller ()
            
            
        end
        
        
        function str = generateMBDynInputString (self)

           str = '';
            
        end
        
    end
    
end