classdef nullDrive < mbdyn.pre.drive
% mbdyn.pre.nullDrive provides a zero-valued drive
%
% Syntax
%
% nd = mbdyn.pre.nullDrive ()
%
% Description
%
% zero valued drive.
%
% mbdyn.pre.nullDrive Methods:
%
%   nullDrive - mbdyn.pre.nullDrive constructor (zero-valued drive)
%   generateMBDynInputString - generates MBDyn input string for the nullDrive
%
    
    methods
        
        function self = nullDrive ()
            % mbdyn.pre.nullDrive constructor (zero-valued drive)
            %
            % Syntax
            %
            % nd = mbdyn.pre.nullDrive ()
            %
            % Description
            %
            % zero valued drive.
            %
            %
            % Output
            %
            %  nd - mbdyn.pre.nullDrive object
            %
            % See Also: 
            %
            
            self.type = 'null';
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for the nullDrive
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (ndrv)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  ndrv - mbdyn.pre.nullDrive object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = self.type;
            
        end
        
    end
    
end