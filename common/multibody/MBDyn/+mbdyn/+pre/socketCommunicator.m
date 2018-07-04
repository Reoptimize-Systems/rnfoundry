classdef socketCommunicator < mbdyn.pre.externalFileCommunicator
% socket based communicator for mbdyn.pre.externalStructuralForce
%
% Description
%
% External file communicator for use with the
% mbdyn.pre.externalStructuralForce element. It specifies the
% use of a socket based communicator to transmit forces and
% kinematics between MBDyn and the external software.
%
% mbdyn.pre.socketCommunicator Methods:
%
%   socketCommunicator - socketCommunicator constructor
%   commInfo - gets communication info for the socket communicator.
%   generateMBDynInputString - generates MBDyn input string
%
%
% See Also: mbdyn.pre.sharedMemoryCommunicator
%
    
    properties (GetAccess = public, SetAccess = private)
        create;
        port;
        path;
        host;
    end
    
    methods
        
        function self = socketCommunicator (varargin)
            % mbdyn.pre.socketCommunicator constructor
            %
            % Syntax
            %
            % sc = socketCommunicator ('Parameter', value)
            %
            % Description
            %
            % External file communicator for use with the
            % mbdyn.pre.externalStructuralForce element. It specifies the
            % use of a socket based communicator to transmit forces and
            % kinematics between MBDyn and the external software.
            %
            % Input
            %
            % Arguments may be supplied as parameter-value pairs. The
            % available options are:
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
            %  'Create' - character vector which can be 'yes' or 'no'. If
            %    'yes', it indicates that MBDyn will create the socket, and
            %    the peer will have to connect to it. Otherwise, when set
            %    to 'no', it indicates that MBDyn will try to connect to an
            %    already existing socket created by the peer. Connecting to
            %    a peer is attempted while reading the input file. Sockets
            %    are created when the input file has been read. MBDyn waits
            %    until all sockets have been connected to by peers before
            %    the simulation is started.
            %
            %  'Host' - MBDyn supports local (unix) sockets, defined using
            %    the Path parameter, and inet sockets, defined using the
            %    Port and Host parameter. When using inet sockets, Host
            %    should be a character vector containing the inet address
            %    of the interface MBDyn will listen on. The use of Host and
            %    Port is mutually exclusive to the use of the Path option.
            %
            %  'Port' - Scalar integer. The Port option is used in
            %    conjuction with the Host option to specify the port inet
            %    sockets will use to communicate on.
            %
            %  'Path' - MBDyn supports local (unix) sockets, defined using
            %    the Path parameter, and inet sockets, defined using the
            %    Port and Host parameter. When using local sockets, Path
            %    should be a character vector containing the file path to
            %    the local socket which MBDyn will listen on.
            %
            % Output
            %
            %  sc - mbdyn.pre.socketCommunicator
            %
            %
            %
            % See Also: 
            %

            options.SleepTime = [];
            options.Coupling = [];
            options.SendAfterPredict = 'yes';
            options.Create = [];
            options.Host = 'localhost';
            options.Port = [];
            options.Path = [];
            
            options = parse_pv_pairs (options, varargin);
            
            
            self = self@mbdyn.pre.externalFileCommunicator ( ...
                        'SleepTime', options.SleepTime, ...
                        'Coupling', options.Coupling, ...
                        'SendAfterPredict', options.SendAfterPredict );
                    
            self.type = 'socket';

            if ~isempty (options.Create)
                self.checkAllowedStringInputs (options.Create, {'yes', 'no'}, true, 'Create');
            end
            
            if ~isempty (options.Path) && ~isempty (options.Port)
                error ('You cannot specify both path and port option for the socket');
            end
            
            if isempty (options.Path) && isempty (options.Port)
                error ('You must specify either Path or Port option for the socket');
            end
            
            if ~isempty (options.Port)
                assert (ischar (options.Host), 'Host must be a character vector');
            end
            
            self.create = options.Create;
            self.path = options.Path;
            self.port = options.Port;
            self.host = options.Host;
            
            if isempty (options.Path)
                self.commMethod = 'inet socket';
            else
                self.commMethod = 'local socket';
            end
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for socket communicator
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (sc)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  sc - mbdyn.pre.socketCommunicator object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = generateMBDynInputString@mbdyn.pre.externalFileCommunicator(self);
            
            if ~isempty (self.create)
                
                str = self.addOutputLine (str, self.commaSepList ('create', self.create), 1, true, 'will mbdyn create socket?');
                
%                 if strcmp (self.Create, 'no')
%                     
%                 end
                
            end
            
            addcomma = any ( ~[ isempty(self.sleepTime);
                                isempty(self.coupling);
                                isempty(self.sendAfterPredict) ] );
                            
            if isempty (self.path)
                % use host and port
                strline = self.commaSepList ('port', self.formatInteger (self.port), 'host', ['"', self.host, '"']); 
                str = self.addOutputLine (str, strline, 1, addcomma);
            else
                % use path
                str = self.addOutputLine (str, self.commaSepList ('path', ['"', self.path, '"']), 1, addcomma);
            end
            
            if ~isempty (self.sleepTime)
                addcomma = any ( ~[ isempty(self.coupling);
                                    isempty(self.sendAfterPredict) ] );
                                
                str = self.addOutputLine (str, self.commaSepList ('sleep time', self.sleepTime), 1, addcomma);
            end
            
            if ~isempty (self.coupling)
                addcomma = ~isempty(self.sendAfterPredict);
                                
                str = self.addOutputLine (str, self.commaSepList ('coupling', self.coupling), 1, addcomma);
            end
            
            if ~isempty (self.sendAfterPredict)
                str = self.addOutputLine (str, self.commaSepList ('send after predict', self.sendAfterPredict), 1, false);
            end
            
        end
        
        function comminfo = commInfo (self)
            % gets communication info for the socket communicator.
            %
            %
            
            comminfo.commMethod = self.commMethod;
            
            if strcmp (self.commMethod, 'local socket')
                comminfo.path = self.path;
            elseif strcmp (self.commMethod, 'inet socket')
                comminfo.host = self.host;
                comminfo.port = self.port;
            else
                error ('unrecognised communication type');
            end

        end
        
    end
    
end