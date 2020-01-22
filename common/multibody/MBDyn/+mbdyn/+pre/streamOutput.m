classdef streamOutput < mbdyn.pre.force
    
    properties (GetAccess = public, SetAccess = protected)
        
        path;
        port;
        host;
        create;
        signal;
        blocking;
        sendFirst;
        abortIfBroken;
        socketType;
        outputSteps;
        echo;
        echoPrecision;
        echoShift;
        streamName;
        
    end
    
    properties (GetAccess = protected, SetAccess = protected)

    end
    
    methods
        
        function self = streamOutput (content, varargin)
            % element to send simulation output to an external program
            %
            % Syntax
            %
            % so = streamOutput (content)
            % so = streamOutput (..., 'Parameter', value)
            %
            % Description
            %
            % element which allows sending simulation data to an external
            % software.
            %
            % Input
            %
            %  content - mbdyn.pre.valuesOutputContent object or an
            %   mbdyn.pre.motionOutputContent object.
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'Create' - true/false flag. If true, it indicates that MBDyn
            %    will create the socket, and the peer will have to connect
            %    to it. Otherwise, when set to false, it indicates that
            %    MBDyn will try to connect to an already existing socket
            %    created by the peer. Connecting to a peer is attempted
            %    while reading the input file. Sockets are created when the
            %    input file has been read. MBDyn waits until all sockets
            %    have been connected to by peers before the simulation is
            %    started.
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
            %  'SendFirst' - true/false flag indicating whether send occurs
            %    before the first time step (by default, data are always
            %    sent)
            %
            %  'AbortIfBroken' - true/false flag indicating whether the
            %    simulation should abort if the connection breaks. If
            %    false, no further data send will occur for the duration of
            %    the simulation (the default).
            %
            %  'SocketType' - character vector indicating the type of
            %    socket to be used. Can be 'tcp' or 'udp'. Default is
            %    'tcp'.
            %
            %  'OutputSteps' - scalar integer indicating that output is to
            %    be sent only every n steps where n is the supplied value
            %    in OutputSteps.
            %
            %  'StreamName' - character vector containing the name of the
            %    RTAI mailbox where the output is written (a unique string
            %    no more than six characters long). It is basically ignored
            %    by the stream output element except when using RTAI.
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
            % Output
            %
            %  so - mbdyn.pre.streamOutput object
            %
            %
            %
            % See Also: mbdyn.mint.MBCNodal
            %
            
            options.Create = [];
            options.Signal = [];
            options.Blocking = [];
            options.SendFirst = [];
            options.AbortIfBroken = [];
            options.SocketType = 'tcp';
            options.OutputSteps = [];
            options.Echo = [];
            options.EchoPrecision = [];
            options.EchoShift = [];
            options.StreamName = sprintf ('%06d', randi (999999)); % only used by RTAI mailbox
            options.Path = [];
            options.Port = [];
            options.Host = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self.type = 'stream output';
            
            assert ( isa (content, mbdyn.pre.valuesOutputContent) ...
                       || isa (content, mbdyn.pre.motionOutputContent), ...
                     'content must be an mbdyn.pre.valuesOutputContent object or an mbdyn.pre.motionOutputContent object' );
            
            if ~isempty (options.Create)
                self.checkLogicalScalar (options.Create, true, 'Create');
            end
            
            if ~isempty (options.Signal)
                self.checkLogicalScalar (options.Signal, true, 'Signal');
            end
            
            if ~isempty (options.Blocking)
                self.checkLogicalScalar (options.Blocking, true, 'Blocking');
            end
            
            if ~isempty (options.SendFirst)
                self.checkLogicalScalar (options.SendFirst, true, 'SendFirst');
            end
            
            if ~isempty (options.AbortIfBroken)
                self.checkLogicalScalar (options.AbortIfBroken, true, 'AbortIfBroken');
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
                    
                    if isempty (options.Create) || strcmpi (options.Create, 'no'), ...
                        error ('If Create is ''no'' or empty, Path cannot be ''auto''');
                    end
                    
                    options.Path = fullfile (tempdir, sprintf ('mbdyn_output_%d.sock', randi (100000) ));
                    
                end
            end
            
            self.create = options.Create;
            self.signal = options.Signal;
            self.blocking = options.Blocking;
            self.sendFirst = options.SendFirst;
            self.abortIfBroken = options.AbortIfBroken;
            self.socketType = options.SocketType;
            self.outputSteps = options.OutputSteps;
            self.path = options.Path;
            self.port = options.Port;
            self.host = options.Host;
            self.echo = options.Echo;
            self.echoPrecision = options.EchoPrecision;
            self.echoShift = options.EchoShift;
            self.streamName = options.StreamName;
            
        end
        
        function str = generateMBDynInputString (self)
            % generates MBDyn input string for streamOutput object
            % 
            % Syntax
            %  
            % str = generateMBDynInputString (so)
            %  
            % Description
            %  
            % generateMBDynInputString is a method shared by all MBDyn
            % components and is called to generate a character vector used
            % to construct an MBDyn input file.
            %  
            % Input
            %  
            %  so - mbdyn.pre.streamOutput object
            %  
            % Output
            %  
            %  str - character vector for insertion into an MBDyn input
            %   file.
            %
            
            str = sprintf ('    %s,', self.type);
            
            if ~isempty (self.create)
                if self.create
                    createstr = 'yes';
                else
                    createstr = 'no'; 
                end
                str = self.addOutputLine (str, self.commaSepList ('create', createstr), 2, true);
            end
            
            if isempty (self.path)
                % use host and port
                strline = self.commaSepList ('port', self.formatInteger (self.port), 'host', ['"', self.host, '"']); 
                str = self.addOutputLine (str, strline, 1, addcomma);
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
            
            if ~isempty (self.sendFirst)
                if self.sendFirst
                    sendFirststr = 'send first';
                else
                    sendFirststr = 'no send first'; 
                end
                str = self.addOutputLine (str, sendFirststr, 2, true);
            end
            
            if ~isempty (self.abortIfBroken)
                if self.abortIfBroken
                    abortIfBrokenstr = 'abort if broken';
                else
                    abortIfBrokenstr = 'do not abort if broken'; 
                end
                str = self.addOutputLine (str, abortIfBrokenstr, 2, true);
            end
            
            if ~isempty (self.outputSteps)
                str = self.addOutputLine (str, self.commaSepList ('output every', self.outputSteps), 2, true);
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
            
            str = self.addOutputLine (str, self.content.generateMBDynInputString (), 2, false);
            
            str = self.addOutputLine (str, ';', 1, false, sprintf ('end %s', self.type));
            
        end
        
    end
    
end