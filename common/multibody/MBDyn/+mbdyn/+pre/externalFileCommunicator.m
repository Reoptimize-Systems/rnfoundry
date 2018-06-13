classdef externalFileCommunicator < mbdyn.pre.base
% base class of the external file communicator classes
%
% Syntax
%
% efc = mbdyn.pre.externalFileCommunicator ()
% efc = mbdyn.pre.externalFileCommunicator ('Parameter', Value)
%
% Description
%
% mbdyn.pre.externalFileCommunicator is a base class for the
% comunicator classes used for e.g. external structural forces.
% COmmunicators define how the passing of data between
% Matlab/Octave and MBDyn is achieved.
%
% mbdyn.pre.externalFileCommunicator Methods:
%
%   externalFileCommunicator - construct an mbdyn.pre.externalFileCommunicator object
%   generateMBDynInputString - generates MBDyn input string for external file communicators
%
%
% See also: mbdyn.pre.socketCommunicator, 
%           mbdyn.pre.sharedMemoryCommunicator
%
    
    properties (GetAccess = public, SetAccess = protected)
        sleepTime;
        coupling;
        precision;
        sendAfterPredict;
        commMethod;
    end
    
    methods
        
        function self = externalFileCommunicator (varargin)
            % construct an mbdyn.pre.externalFileCommunicator object
            %
            %
            % Syntax
            %
            % efc = mbdyn.pre.externalFileCommunicator ()
            % efc = mbdyn.pre.externalFileCommunicator ('Parameter', Value)
            %
            % Description
            %
            % mbdyn.pre.externalFileCommunicator is a base class for the
            % comunicator classes used for e.g. external structural forces.
            % COmmunicators define how the passing of data between
            % Matlab/Octave and MBDyn is achieved.
            %
            % Input
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'SleepTime' - determines how long MBDyn is supposed to sleep
            %    while waiting for a new input file to appear or for an old
            %    output file to disappear.
            %
            %  'Coupling' - character vector determining the type of
            %    coupling with the external software. Can be one of
            %    'staggered', 'tight' or 'loose'. For details on the
            %    implications of each of these, see the section on External
            %    Force in the official MBDyn manual.
            %
            %  'SendAfterPredict' - character vector which can be 'yes' or
            %    'no'. This indicates whether MBDyn must send the predicted
            %    motion or not when playing a tight coupling loop.
            %
            %  'Precision' - only relevant for some communicators. Defines
            %    precision with which data will be passed.
            %
            % Output
            %
            %  efc - mbdyn.pre.externalFileCommunicator
            %
            %
            %
            % See Also: mbdyn.pre.socketCommunicator, 
            %           mbdyn.pre.sharedMemoryCommunicator
            %

            options.SleepTime = [];
            options.Precision = [];
            options.Coupling = 'loose';
            options.SendAfterPredict = 'yes';
            
            options = parse_pv_pairs (options, varargin);
            
            if ~isempty (options.Coupling)
                if ischar (options.Coupling)
                    
                    self.checkAllowedStringInputs ( options.Coupling, ...
                        {'staggared', 'loose', 'tight'}, ...
                        true, 'Coupling' );
                    
                elseif ~self.checkScalarInteger (options.Coupling, false)
                    
                    error ('Coupling must be a string or integer number of steps');
                    
                end
            end
            
            if ~( ( isnumeric (options.SleepTime) ...
                        && isscalar (options.SleepTime) ...
                        && options.SleepTime >= 0 ) ...
                    || isempty (options.SleepTime) )
                
                error ('SleepTime must be a scalar numeric value >= 0');
                
            end
            
            if ~( ( self.checkScalarInteger (options.Precision, false) ) ...
                  || isempty (options.Precision) )
                
                error ('Precision must be an integer');
                
            end
            
            if ~isempty (options.SendAfterPredict)
                self.checkAllowedStringInputs (options.SendAfterPredict, {'yes', 'no'}, true, 'SendAfterPredict');
            end
            
            self.sleepTime = options.SleepTime;
            self.coupling = options.Coupling;
            self.precision = options.Precision;
            self.sendAfterPredict = options.SendAfterPredict;
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for external file communicators
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (efc)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  efc - mbdyn.pre.externalFileCommunicator object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = sprintf ('%s,', self.type);
            
        end
        
%         function comminfo = commInfo (self)
%             % gets communication info for the socket communicator.
%             %
%             %
%             
%             comminfo.commMethod = self.commMethod;
%             
%             if strcmp (self.commMethod, 'local socket')
%                 comminfo.path = self.path;
%             elseif strcmp (self.commMethod, 'inet socket')
%                 comminfo.host = self.host;
%                 comminfo.port = self.port;
%             else
%                 error ('unrecognised communication type');
%             end
% 
%         end
        
    end
    
    
end