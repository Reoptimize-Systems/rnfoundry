classdef scalarFunction < mbdyn.pre.base
    
    properties (GetAccess = public, SetAccess = protected)
        name;
        fcnType;
    end
    
    methods
        
        function self = scalarFunction (name, fcnType)
            % base class for scalar function objects
            
            assert (ischar (name), 'name should be a char array');
            assert (ischar (fcnType), 'fcnType should be a char array');
            
            self.type = 'scalar function';
            self.name = name;
            self.fcnType = fcnType;
            
        end
        
%         function str = generateOutputString (self)
%             
%            
%         end
        
    end
    
end