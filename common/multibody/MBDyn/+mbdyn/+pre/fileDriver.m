classdef fileDriver < mbdyn.pre.driver
% base class for all MBDyn file drivers
%
% Syntax
%
% fd = mbdyn.pre.fileDriver ('Parameter', Value)
%
% Description
%
% mbdyn.pre.fileDriver is the base class for all other fileDriver types in
% the toolbox. It contains methods and properties common to all
% fileDrivers. It is not intended to be used directly by ordinary users.
%
% mbdyn.pre.fileDriver Methods:
%
%   fileDriver - mbdyn.pre.fileDriver constructor
%   
%
    
    properties (GetAccess = public, SetAccess = public)
       
        subType; % subType of the file driver
        
    end
    
    properties (GetAccess = public, SetAccess = protected)
       

        
    end
    
    properties (GetAccess = protected, SetAccess = protected)
       

        
    end
    
    methods
        
        function self = fileDriver ()
            % mbdyn.pre.fileDriver constructor
            %
            % Syntax
            %
            % fd = mbdyn.pre.fileDriver ('Parameter', Value)
            %
            % Description
            %
            % mbdyn.pre.fileDriver is the base class for all other fileDriver
            % types in the toolbox. It contains methods and properties
            % common to all fileDrivers such as bodies, joints and forces etc.
            % It is not intended to be used directly by ordinary users.
            %
            % Input
            %
            % Arguments may be supplied as parameter-value pairs. The
            % available options are:
            %
            %
            % Output
            %
            %  fd - mbdyn.pre.fileDriver object
            %
            
            self = self@mbdyn.pre.driver ();
            
            self.type = 'file';
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for a file driver
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (fd)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  fd - mbdyn.pre.fileDriver object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = self.addOutputLine ('' , sprintf ('%s : %d, %s,', self.type, self.label, self.subType), 1, false, sprintf ('file %s driver', self.subType));
            
            % delete newline character and space from start
            str(1) = [];
            
        end
        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = struct ();
            
            nopass_list = {};
            
        end
        
    end
    
end