classdef stringDrive < mbdyn.pre.drive
    
    properties (GetAccess = public, SetAccess = protected)
        
        string;
        labelReplacementObjects;
        
    end
    
    methods
        
        function self = stringDrive (string, varargin)
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
            
            [options, nopass_list] = mbdyn.pre.stringDrive.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            pvpairs = mbdyn.pre.base.passThruPVPairs (options, nopass_list);
            
            self = self@mbdyn.pre.drive (pvpairs{:});

            if ~ischar (string)
                error ('''string'' must be a char array')
            end
            
            self.string = string;
            self.type = 'string';
            self.labelReplacementObjects = options.LabelRepObjects;
                
        end
        
        function str = generateMBDynInputString (self)
            
            drivestr = self.string;
            
            if ~isempty (self.labelReplacementObjects)
%                 uids = regexp (self.string, 'UID:([0-9]+)', 'tokens', 'forceCellOutput');
%                 uids = uids{1};
                for objind = 1:numel (self.labelReplacementObjects)
                    drivestr = strrep ( drivestr, ...
                                        ['UID:', int2str(self.labelReplacementObjects{objind}.uid)], ...
                                        int2str(self.labelReplacementObjects{objind}.label) );
                end
                
            end
            
            str = self.commaSepList (self.type, ['"', drivestr, '"']);
            
        end
        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options = mbdyn.pre.drive.defaultConstructorOptions ();
            
            parentfnames = fieldnames (options);
            
            options.LabelRepObjects = {};
            
            allfnames = fieldnames (options);
            
            nopass_list = setdiff (allfnames, parentfnames, 'stable');
            
        end
        
    end
    
end