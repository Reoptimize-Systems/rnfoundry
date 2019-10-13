classdef fileDrive < mbdyn.pre.drive
    
    properties (GetAccess = public, SetAccess = private)
        fileDriver;
        columnNumber;
        amplitude;
    end
    
    methods
        
        function self = fileDrive (file_driver, varargin)
            % gets file driver input
            %
            % Syntax
            %
            % ld = fileDrive (const_coef, slope_coef)
            %
            % Description
            %
            % This drive is attached to a file driver object (defined in
            % the 'Drivers' section of an MBDyn input file). It is used to
            % provide the data from an external program to elements which
            % can use a driveCaller.
            %
            % Input
            %
            %  file_driver - mbdyn.pre.fileDriver object
            %
            %
            % Output
            %
            %  ld - mbdyn.pre.fileDrive object
            %
            % See Also: 
            %
            
            [options, ~] = mbdyn.pre.fileDrive.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            

            assert ( isa (file_driver, 'mbdyn.pre.fileDriver'), ...
                     'file_driver must be an mbdyn.pre.fileDriver object' );
                 
            if ~isempty (self.columnNumber)
                self.checkScalarInteger (options.ColumnNumber, true, 'ColumnNumber');
            end
            
            if ~isempty (self.amplitude)
                self.checkNumericScalar (options.Amplitude, true, 'Amplitude');
            end
            
            self.fileDriver = file_driver;
            self.columnNumber = options.ColumnNumber;
            self.amplitude = options.Amplitude;
            
            self.type = 'file';
            
        end
        
        function str = generateMBDynInputString (self)
            
            str = self.commaSepList (self.type, mbdyn.pre.base.formatInteger (self.fileDriver.label));
            
            if ~isempty (self.columnNumber)
                str = sprintf ('%s, %s', str, mbdyn.pre.base.formatInteger (self.columnNumber) );
            end
            
            if ~isempty (self.amplitude)
                str = self.commaSepList (str, self.amplitude);
            end
            
        end
        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options.ColumnNumber = [];
            options.Amplitude = [];
            
            nopass_list = {};
            
        end
        
    end
    
end