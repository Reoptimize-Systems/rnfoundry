classdef directDrive < mbdyn.pre.drive
    
    properties (GetAccess = public, SetAccess = private)

    end
    
    methods
        
        function self = directDrive ()
            % Constructor for direct drive
            %
            % Syntax
            %
            % dd = mbdyn.pre.directDrive ()
            %
            % Description
            %
            % Transparently returns the input value; the arglist is empty.
            % It is useful in conjunction with those drive callers that
            % require their output to be fed into another drive caller,
            % like the dof, node and element drive callers, when the output
            % needs to be used as is.
            %
            %
            % Output
            %
            %  dd - mbdyn.pre.directDrive object
            %
            %
            %
            % See Also: 
            %
            
            self.type = 'direct';
            
        end
        
        function str = generateOutputString (self)
            
            str = self.type;
            
        end
        
    end
    
end