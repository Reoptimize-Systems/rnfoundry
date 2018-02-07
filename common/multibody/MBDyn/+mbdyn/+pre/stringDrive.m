classdef stringDrive < mbdyn.pre.drive
    
    properties (GetAccess = public, SetAccess = protected)
        
        string;
        
    end
    
    methods
        
        function self = stringDrive (string)
            % drive which returns the value of a mathematical expression
            %
            % Syntax
            %
            % sd = stringDrive (string)
            %
            % Description
            %
            % stringDrive is a drive which represents a mathematical
            % expression which is evaluated at each time step of the
            % simulation. 
            %
            % Input
            %
            %  string - expression parsed by the math parser every time the
            %    drive is invoked. Two special variable may be used in the
            %    expression 'Time' which is the current simulation time,
            %    and 'Var'. The value of 'Var' is set by the caller of the
            %    string drive
            %
            % Output
            %
            % sd - mbdyn.pre.stringDrive object
            %
            % See Also: 
            %

            if ~ischar (string)
                error ('''string'' must be a char array')
            end
            
            self.string = string;
                
        end
        
        function str = generateOutputString (self)
            
            str = self.commaSepList (self.type, ['"', self.string, '"']);
            
        end
        
    end
    
end