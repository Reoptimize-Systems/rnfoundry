classdef multDrive < mbdyn.pre.drive
    
    properties (GetAccess = public, SetAccess = private)
        drive1;
        drive2;
    end
    
    methods
        
        function self = multDrive (drive1, drive2)
            % Constructor for mult drive
            %
            % Syntax
            %
            % md = mbdyn.pre.multDrive (drive1, drive2)
            %
            % Description
            %
            % multDrive multiplies the value of two subordinate drives.
            %
            % Input
            %
            %  drive1 - 
            %
            %  drive2 - 
            %
            % Output
            %
            %  md - mbdyn.pre.multDrive object
            %
            %
            %
            % See Also: 
            %
            
            assert (isa (drive1, 'mbdyn.pre.drive'), ...
                'drive1 must be an mbdyn.pre.drive' );
            
            assert (isa (drive2, 'mbdyn.pre.drive'), ...
                'drive2 must be an mbdyn.pre.drive' );
            
            self.drive1 = drive1;
            self.drive2 = drive2;
            self.type = 'mult';
            
        end
        
        function str = generateOutputString (self)
            
            str = [ self.type, ',' ];
            
            str = self.addOutputLine ( str, ...
                self.drive1.generateOutputString (), ...
                                       1, ...
                                       true );
                                   
            str = self.addOutputLine ( str, ...
                self.drive2.generateOutputString (), ...
                                       1, ...
                                       false );
            
        end
        
    end
    
end