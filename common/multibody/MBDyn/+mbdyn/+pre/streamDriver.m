classdef streamDriver < mbdyn.pre.fileDriver
% base class for all MBDyn file drivers
%
% Syntax
%
% fd = mbdyn.pre.streamDriver ('Parameter', Value)
%
% Description
%
% mbdyn.pre.streamDriver is the base class for all other streamDriver types in
% the toolbox. It contains methods and properties common to all
% fileDrivers. It is not intended to be used directly by ordinary users.
%
% mbdyn.pre.streamDriver Methods:
%
%   streamDriver - mbdyn.pre.streamDriver constructor
%   
%
    
    
    properties (GetAccess = public, SetAccess = protected)
       
        path;
        port;
        host;
        create;
        signal;
        blocking;
        receiveFirst;
        socketType;
        inputSteps;
        echo;
        echoPrecision;
        echoShift;
        streamName;
        timeout;
        commMethod;
        nValues;
        initialValues;
        
    end
    
    properties (GetAccess = protected, SetAccess = protected)
       
        socket;
        outputStream;
        dataOutputStream;
        socketReady;
        socketNBytes;
        
    end
    
    methods
        
        function self = streamDriver (create, nvalues, varargin)
            % element to send simulation output to an external program
            %
            % Syntax
            %
            % so = streamDriver (content)
            % so = streamDriver (..., 'Parameter', value)
            %
            % Description
            %
            % force element which allows communication with an external
            % software that computes forces applied to a pool of nodes and
            % may depend on the kinematics of those nodes.
            %
            % Input
            %
            %  content - mbdyn.pre.valuesOutputContent object or an
            %   mbdyn.pre.motionOutputContent object.
            %
            %  create - true/false flag. If true, it indicates that MBDyn
            %   will create the socket, and the peer will have to connect
            %   to it. Otherwise, when set to false, it indicates that
            %   MBDyn will try to connect to an already existing socket
            %   created by the peer. Connecting to a peer is attempted
            %   while reading the input file. Sockets are created when the
            %   input file has been read. MBDyn waits until all sockets
            %   have been connected to by peers before the simulation is
            %   started.
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'StreamName' - character vector containing the name of the
            %    RTAI mailbox where the output is written (a unique string
            %    no more than six characters long). It is basically ignored
            %    by the stream output element except when using RTAI. If
            %    not supplied a random integer between 1 and 999999 is
            %    chosen and converted to a string.
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
            %    should be a character vector. It can either contain the
            %    file path to the local socket which MBDyn will listen on,
            %    or can be the string 'auto'. If it is 'auto' a socket file
            %    name will be generated with a random integer in the name
            %    in the temp directory with the pattern:
            %
            %    /tmp/mbdyn_<random_integer>.sock
            %
            %    e.g. /tmp/mbdyn_172727.sock
            %
            %    The integer will be generated using the randi function and
            %    will be between 1 and 100000. The result is stored in the
            %    path proerty for later use. If the path is set to 'auto'
            %    and 'Create' is 'no' an error will be thrown.
            %
            %  'Signal' - true/false flag. If true SIGPIPE will be raised
            %    when sending through a socket when the other end is
            %    broken, if false SIGPIPE will not be raised (by
            %    default, SIGPIPE is raised).
            %
            %  'Blocking' - true/false flag indicating whether operations
            %    on the socket block. Default is operation do block.
            %
            %  'ReceiveFirst' - true/false flag indicating whether receive
            %    occurs before the initial time step (by default, data are
            %    expected to be received).
            %
            %  'SocketType' - character vector indicating the type of
            %    socket to be used. Can be 'tcp' or 'udp'. Default is
            %    'tcp'.
            %
            %  'InputSteps' - scalar integer indicating that input is to
            %    be received only every n steps where n is the supplied value
            %    in InputSteps.
            %
            %  'Echo' - optional character vector containing a path to a
            %    file where the output of the reference configuration will
            %    be exported. If the file exists, it is overwritten, if
            %    allowed by file system permissions. The format is that of
            %    the communicator in stream form. If empty nothing will be
            %    exported. Default is empty.
            %
            %  'EchoPrecision' - optional parameter determining the number
            %    of digits used in Echo output file.
            %
            %  'Timeout' - 
            %
            %
            % Output
            %
            %  so - mbdyn.pre.streamDriver object
            %
            %
            %
            % See Also: mbdyn.mint.MBCNodal
            %
            
            [options, ~] = mbdyn.pre.streamDriver.defaultConstructorOptions ();
            
            options = parse_pv_pairs (options, varargin);
            
            self = self@mbdyn.pre.fileDriver ();

            self.subType = 'stream';
            
            self.checkLogicalScalar (create, true, 'create');
            
            self.checkScalarInteger (nvalues, true, 'nvalues');
            
            if ~isempty (options.Signal)
                self.checkLogicalScalar (options.Signal, true, 'Signal');
            end
            
            if ~isempty (options.Blocking)
                self.checkLogicalScalar (options.Blocking, true, 'Blocking');
            end
            
            if ~isempty (options.ReceiveFirst)
                self.checkLogicalScalar (options.ReceiveFirst, true, 'ReceiveFirst');
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
            
            if isempty (options.Path)
                self.commMethod = 'inet socket';
            else
                self.commMethod = 'local socket';
                
                if strcmpi (options.Path, 'auto')
                    
                    if isempty (create) || strcmpi (create, 'no'), ...
                        error ('If Create is ''no'' or empty, Path cannot be ''auto''');
                    end
                    
                    options.Path = fullfile (tempdir, sprintf ('mbdyn_output_%d.sock', randi (100000) ));
                    
                end
            end
            
            if ~isempty (options.InitialValues)
                
                assert ( isvector (options.InitialValues) ...
                           && numel (options.InitialValues) == nvalues, ...
                         'InitialValues must be a real numeric vector of length nvalues' );
                     
                if iscell (options.InitialValues)
                    for ind = 1:numel (options.InitialValues)
                        self.checkNumericScalar ( options.InitialValues{ind}, true, sprintf ('InitialValues{%d}', ind));
                    end
                elseif isnumeric (options.InitialValues)
                    % convert it to a cell array
                    options.InitialValues = mat2cell (options.InitialValues(:), ones (1, nvalues));
                else
                    error ('InitialValues must be a real numeric vector of length nvalues or a cell array of length nvalues, each cell containing scalar values' );
                end
                
            end
            
            if ~isempty (options.Timeout)
                self.checkNumericScalar (options.Timeout, true, '1Timeout');
                assert (options.Timeout >= 0, 'Timeout must be greatr than or equal to zero')
            end
            
            self.create = create;
            self.nValues = nvalues;
            self.signal = options.Signal;
            self.blocking = options.Blocking;
            self.receiveFirst = options.ReceiveFirst;
            self.socketType = options.SocketType;
            self.inputSteps = options.InputSteps;
            self.path = options.Path;
            self.port = options.Port;
            self.host = options.Host;
            self.echo = options.Echo;
            self.echoPrecision = options.EchoPrecision;
            self.echoShift = options.EchoShift;
            self.streamName = options.StreamName;
            self.socketReady = false;
            self.initialValues = options.InitialValues;
            self.timeout = options.Timeout;
            
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
            %  fd - mbdyn.pre.streamDriver object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = generateMBDynInputString@mbdyn.pre.fileDriver (self);
            
            str = self.addOutputLine (str, self.commaSepList ('name', sprintf ('"%s"', self.streamName)), 2, true);
            
            if self.create
                createstr = 'yes';
            else
                createstr = 'no'; 
            end
            str = self.addOutputLine (str, self.commaSepList ('create', createstr), 2, true);
            
            if isempty (self.path)
                % use host and port
                strline = self.commaSepList ('port', self.formatInteger (self.port), 'host', ['"', self.host, '"']); 
                str = self.addOutputLine (str, strline, 2, true);
            else
                % use path
                str = self.addOutputLine (str, self.commaSepList ('local', ['"', self.path, '"']), 2, true);
            end
            
            if ~isempty (self.socketType)
                str = self.addOutputLine (str, self.commaSepList ('socket type', self.socketType), 2, true);
            end
            
            if ~isempty (self.signal)
                if self.signal
                    signalstr = 'signal';
                else
                    signalstr = 'no signal'; 
                end
                str = self.addOutputLine (str, signalstr, 2, true);
            end

            if ~isempty (self.blocking)
                if self.blocking
                    blockingstr = 'blocking';
                else
                    blockingstr = 'non blocking'; 
                end
                str = self.addOutputLine (str, blockingstr, 2, true);
            end

            if ~isempty (self.inputSteps)
                str = self.addOutputLine (str, self.commaSepList ('input every', self.inputSteps), 2, true);
            end
            
            if ~isempty (self.receiveFirst)
                if self.receiveFirst
                    rcvFirststr = 'yes';
                else
                    rcvFirststr = 'no'; 
                end
                str = self.addOutputLine (str, self.commaSepList ('receive first', rcvFirststr), 2, true);
            end
            
            if ~isempty (self.timeout)
                str = self.addOutputLine (str, self.commaSepList ('timeout', self.timeout), 2, true, 'echo output to file');
            end

            if ~isempty (self.echo)
                str = self.addOutputLine (str, self.commaSepList ('echo', ['"', self.echo, '"']), 2, true, 'echo output to file');
            end
            
            if ~isempty (self.echoPrecision)
                str = self.addOutputLine (str, self.commaSepList ('precision', self.echoPrecision), 2, true, 'echo output to file');
            end
            
            if ~isempty (self.echoShift)
                str = self.addOutputLine (str, self.commaSepList ('shift', self.echoShift), 2, true, 'echo output to file');
            end
            
            addcomma = ~isempty (self.initialValues);
            
            str = self.addOutputLine (str, sprintf ('%d', self.nValues), 2, addcomma);
            
            if ~isempty (self.initialValues)
                 str = self.addOutputLine (str, self.commaSepList ('initial values', self.initialValues{:}), 2, false, 'initial values');
            end
            
            % TODO: stream input content modifiers
%             str = self.addOutputLine (str, self.content.generateMBDynInputString (), 2, false);
            
            str = self.addOutputLine (str, ';', 1, false, sprintf ('end %s %s', self.subType, self.type));
            
        end
        
        function status = start (self)
            
            if ~isempty (self.path)
                error ('Unix (local) sockets are not yet supported');
            end
            
            status = -1;
            
            % create input socket (2 double: x, x_prime)
            if self.create
                % if self.create is true, it means MBDyn will create the
                % socket, and we will connect to it here
                try
                    self.socket = comms.jtcp('request', self.host, self.port, 'timeout', 60e3, 'serialize', false);
                catch
                    return;
                end
                
                status = 0;
            else
                % self.create is false we need to make the socket, MBDyn
                % will not create it
                self.socket = comms.jtcp('accept', self.port, 'timeout', 60e3, 'serialize', false);
            end
            
            self.socketNBytes = self.nValues * 8;
            
%             self.socket = java.net.Socket(self.host, self.port);
            
%             self.ouputStream = self.socket.getOutputStream ();
            
%             self.dataOutputStream = java.io.DataOutputStream(self.ouputStream);
            
            self.socketReady = true;
            
        end
        
        function sendValues (self, values)
            
            assert (self.socketReady, 'Cannot send values socket is not ready yet. Have you called startSocket?');
            
            assert ( numel (values) == self.nValues, ...
                     'numel(values) should be %d, but was %d', self.nValues, numel (values) );
            
            for ind = 1:numel (values)
                comms.jtcp('write', self.socket, typecast (values{ind}, 'int8'));
            end
            
%             self.dataOutputStream.writeBytes(typecast(values, 'unint8'));
            
%             self.dataOutputStream.flush;
            
        end
        
        function stop (self)
            
            self.socket.close ();
            
            self.socketReady = false;
            
        end
        
    end
    
    methods (Static)
        
        function [options, nopass_list] = defaultConstructorOptions ()
            
            options.Signal = [];
            options.Blocking = [];
            options.ReceiveFirst = [];
            options.AbortIfBroken = [];
            options.SocketType = 'tcp';
            options.InputSteps = [];
            options.Echo = [];
            options.EchoPrecision = [];
            options.EchoShift = [];
            options.StreamName = sprintf ('%06d', randi (999999)); % only used by RTAI mailbox
            options.Path = [];
            options.Port = [];
            options.Host = [];
            options.InitialValues = [];
            options.Timeout = [];
            
            nopass_list = {};
            
        end
        
    end
    
end