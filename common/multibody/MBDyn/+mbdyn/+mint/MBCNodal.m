classdef MBCNodal < mbdyn.mint.cppinterface
    % MBCNodal: interface for MBDyn multibody dynamics library using sockets. 
    %
    % Description
    %
    % 

    properties (SetAccess = private, GetAccess = public)
        
        NNodes; % Number of structural external force nodes
        structuralNodes; % cell array of structural nodes from the structural external force object (if supplied)
        
        useLabels; % Flag indicating whether to make node label numbers available
        useAccelerations; % Flag indicating whether accelerations are to be made available
        useRefNode; % Flag indicating whether there is a reference node
        
        useMoments; % Flag indicating whether moments are to be used
        dataAndNext;
        
        
        host;
        port;
        path;
        sharedMemoryName;
        commMethod;
        MBDynInputFile;
        MBDynExecutable;
        MBDynStartWaitTime;
        MBDynOutputFile;
        
    end
    
    properties (SetAccess = private, GetAccess = private)
        
        % needMoments -Internal flag used to determine if we need to apply
        % moments before advancing
        needMoments;
        
        % needForces -Internal flag used to determine if we need to apply
        % forces before advancing 
        needForces;
        
        mbsys;
        outputPrefix;
        nodeOrientiationType;
        setNodePositions;
        
    end
    
	% enum not yet implemented in Octave
%     enumeration
%         NOT_READY,
% 		INITIALIZED,
% 		SHARED_MEMORY_READY,
% 		READY,
% 		FINISHED
%     end
    
    methods
        
        % Constructor
        function self = MBCNodal (varargin)
            % Creates an MBCNodal object
            %
            % Syntax
            %
            % mb = MBCNodal ('Parameter', value)
            %
            % Input
            % 
            % Arguments are be supplied using parameter-value pair syntax.
            % Which arguments are required and which are optional depends
            % on how the system is defined.
            %
            % 'MBDynPreProc' - use this option to supply an
            %   mbdyn.pre.system object from which the communication data
            %   will be determined. If this options is used, the optional
            %   arguments 'NNodes', 'NodeOrientationType', 'UseLabels',
            %   'UseAccelerations', 'UseRefNode' and
            %   'RefNodeOrientationType' will be ignored as they will be
            %   set through examination of the external structual force
            %   element in the MBDyn system description. Note that
            %   'UseMoments' will still be used, this is not set in the
            %   input file.
            %
            % 'MBDynInputFile' - MBDyn input file name. This is the file
            %   containing the problem description for MBDyn. This can
            %   either be a pre-existing file, or can be created by
            %   MBCNodal from a mbdyn.pre.system supplied using the
            %   'MBDynPreProc' option. If no mbdyn.pre.system object is
            %   supplied, this must be a preexisting file.
            %
            % 'CreateInputFile' - if a mbdyn.pre.system object is supplied
            %   using the MBDynPreProc option, the input mbd file for MBDyn
            %   can be created from the system. This options is a logical
            %   flag determining whether the file should be generated. If
            %   no file path is supplied via the 'MBDynInputFile' option, a
            %   random temporary file name will be generated. The name of
            %   the resulting file will be stored in the MBDynInputFile
            %   class property. Default is true.
            %
            % 'OverwriteInputFile' - if a mbdyn.pre.system object is
            %   supplied, and 'CreateInputFile' is true, this flag
            %   determines whether to overwrite any existing file in that
            %   location. Default is true. 
            %
            % 'OutputPrefix' - used to specify the output path prefix. This
            %   is the name of output files, but without their file
            %   extension. e.g. /home/jbloggs/my_mbdyn_sim will create the
            %   files:
            %
            %   /home/jbloggs/my_mbdyn_sim.frc
            %   /home/jbloggs/my_mbdyn_sim.ine
            %   /home/jbloggs/my_mbdyn_sim.out
            %   /home/jbloggs/my_mbdyn_sim.mov
            %   /home/jbloggs/my_mbdyn_sim.jnt
            %   /home/jbloggs/my_mbdyn_sim.log
            %
            %   A windows example might look like
            %   C:\Users\IEUser\Documents\my_mbdyn_sim 
            %   producing the files:
            %
            %   C:\Users\IEUser\Documents\my_mbdyn_sim.frc
            %   C:\Users\IEUser\Documents\my_mbdyn_sim.ine
            %   C:\Users\IEUser\Documents\my_mbdyn_sim.out
            %   C:\Users\IEUser\Documents\my_mbdyn_sim.mov
            %   C:\Users\IEUser\Documents\my_mbdyn_sim.jnt
            %   C:\Users\IEUser\Documents\my_mbdyn_sim.log
            %
            % 'MBDynExecutable' - used to specify the full path to the
            %   mbdyn executeable file. If not supplied, MBCNodal will look
            %   in various standard location where it might be found. An
            %   error will be thrown if MBDyn cannot be located. Default is
            %   an empty string, meaning MBCNodal will search for MBDyn.
            %
            % 'UseMoments' - optional logical flag indicating whether data
            %   on moments is to be exchanged. Default is false.
            %
            % 'DataAndNext' - optional logical flag indicating when data
            %   will be exchanged. Default is true.
            %
            % The following options are used when the MBDynPreProc options
            % is not used to supply a mbdyn.pre.system object.
            %
            % 'NNodes' - scalar number of nodes in the problem, should
            %   match up with problem file. This MUST be supplied if you
            %   are not using a reference node.
            %
            % 'NodeOrientationType' - string determining the format in
            %   which node orientations will be returned by MBDyn. Can be
            %   one of 'none', 'theta' , 'mat' and 'euler123'. Default is
            %   'mat'.
            %
            % 'UseLabels' - true/false flag determining whether node labels
            %   will be used and made available. Defautl is true.
            %
            % 'UseAccelerations' - true/false flag determining whether
            %   accelerations will be returned by MBDyn. Default is false.
            %
            % 'UseRefNode' - true/false flag indicating whether a reference
            %   node is to be used. If so, you should provide a value for
            %   the 'RefNodeOrientationType' option. Default is false.
            %
            % 'RefNodeOrientationType' - string denoting the reference node
            %   orientation matrix type to be used. This determines the
            %   format in which node orientations will be returned by
            %   MBDyn. Can be one of: 'none', 'theta' , 'mat' and
            %   'euler123'. If UseRefNode is false, this will be ignored.
            %   If UseRefNode is true and 'RefNodeOrientationType' is
            %   'none', the value supplied in the 'NodeOrientationType'
            %   option will be used. Default is 'none'.
            %
            % 
            
            options.MBDynInputFile = '';
            options.CreateInputFile = true;
            options.OverwriteInputFile = true;
            options.OutputPrefix = '';
            options.MBDynExecutable = '';
            options.UseMoments = false;
            options.DataAndNext = true;
            
            % mbdyn.pre.system is supplied via this
            options.MBDynPreProc = [];
            options.SetStructuralNodePositions = true;
            
            % the following only used if no mbdyn.pre.system is supplied
            options.CommMethod = '';
            options.Host = '';
            options.Port = '';
            options.Path = '';
            options.SharedMemoryName = '';
            options.UseRefNode = false; 
            options.RefNodeOrientationType = 'none'; % refnode rotation type 'none'. Will be ignored if UseRefNode is false
            options.NNodes = 0; % number of nodes, should match up with problem file
            options.NodeOrientationType = 'orientation matrix';
            options.UseLabels = true;
            options.UseAccelerations = false; % don't handle accelerations by default
            
            options = parse_pv_pairs (options, varargin);
            
            assert (ischar (options.MBDynInputFile), 'MBDynInputFile must be a char array');
            mbdyn.pre.base.checkLogicalScalar (options.CreateInputFile, true, 'CreateInputFile');
            mbdyn.pre.base.checkLogicalScalar (options.OverwriteInputFile, true, 'OverwriteInputFile');
            assert (ischar (options.OutputPrefix), 'OutputPrefix must be a char array');
            assert (ischar (options.MBDynExecutable), 'MBDynExecutable must be a char array');
            mbdyn.pre.base.checkLogicalScalar (options.UseMoments, true, 'UseMoments');
            mbdyn.pre.base.checkLogicalScalar (options.DataAndNext, true, 'DataAndNext');
            
            if isa (options.MBDynPreProc, 'mbdyn.pre.system')
                
                comminfo = options.MBDynPreProc.externalStructuralCommInfo ();
                
                if strcmp (comminfo.commMethod, 'inet socket') ...
                    || strcmp (comminfo.commMethod, 'local socket')
                        % initialise the cppinterface parent class by passing the
                        % mexfunction to the superclass constructor
                        mexfcn = @mbdyn.mint.mexMBCNodal;
                elseif strcmp (comminfo.commMethod, 'shared memory')
                        mexfcn = @mbdyn.mint.mexMBCNodalSharedMem;
                else
                    error ('unrecognised communication method');
                end

            else
                
                switch options.CommMethod
                    
                    case 'inet socket'
                        
                        mexfcn = @mbdyn.mint.mexMBCNodal;
                        
                    case 'local socket'
                        
                        mexfcn = @mbdyn.mint.mexMBCNodal;
                        
                    case 'shared memory'
                        
                        mexfcn = @mbdyn.mint.mexMBCNodalSharedMem;
                        
                    otherwise
                        
                end
                
            end
            
            self = self@mbdyn.mint.cppinterface(mexfcn);
            
            if ~isempty (options.MBDynPreProc)
                
                % use the information in the MBDynPreProc class to set up
                % the sim communication
                if isa (options.MBDynPreProc, 'mbdyn.pre.system')
                    
                    self.mbsys = options.MBDynPreProc;
                    
                    comminfo = self.mbsys.externalStructuralCommInfo ();
                    
                    self.commMethod = comminfo.commMethod;
                    
                    if strcmp (comminfo.commMethod, 'inet socket')
                        self.host = comminfo.host;
                        self.port = comminfo.port;
                    elseif strcmp (comminfo.commMethod, 'local socket')
                        self.path = comminfo.path;
                    elseif strcmp (comminfo.commMethod, 'shared memory')
                        self.sharedMemoryName = comminfo.shared_mem_name;
                    end
                    
                else
                    error ('MBDynPreProc is not empty or an mbdyn.pre.system object');
                end
                
                mbdyn.pre.base.checkLogicalScalar ( options.SetStructuralNodePositions, ...
                                                    true, 'SetStructuralNodePositions' );
                
                if ~isempty (options.MBDynInputFile)
                    self.MBDynInputFile = options.MBDynInputFile;
                else
                    [pathstr, name] = fileparts (tempname ());
                    self.MBDynInputFile = fullfile (pathstr, ['mbdyn_input_file_', name, '.mbd']);
                end
                
                if exist (self.MBDynInputFile, 'file') == 2
                    % file exists
                    if options.OverwriteInputFile
                        self.mbsys.generateMBDynInputFile (self.MBDynInputFile);
                    end
                else
                    if options.CreateInputFile
                        self.mbsys.generateMBDynInputFile (self.MBDynInputFile);
                    end
                end
                
                self.useMoments = options.UseMoments;
                
                extfinfo = self.mbsys.externalStructuralInfo ();
                
                % copy some of the chosen options over to the class properties
                self.useLabels = extfinfo.UseLabels;
                self.useAccelerations = extfinfo.UseAccelerations;
                self.useRefNode = extfinfo.UseRefNode;
                self.nodeOrientiationType = extfinfo.NodeOrientationType;
                
                self.NNodes = extfinfo.NNodes;
                self.structuralNodes = extfinfo.Nodes;
                self.setNodePositions = options.SetStructuralNodePositions;
                
            else
                % No mbsys object supplied
                mbdyn.pre.base.checkLogicalScalar (options.UseAccelerations, true, 'UseAccelerations');
                mbdyn.pre.base.checkLogicalScalar (options.UseLabels, true, 'UseLabels');
                mbdyn.pre.base.checkLogicalScalar (options.UseRefNode, true, 'UseRefNode');
                
                
                if ~isempty (options.MBDynInputFile)
                    if exist (options.MBDynInputFile, 'file') == 2
                        self.MBDynInputFile = options.MBDynInputFile;
                    else
                        error ('MBDynInputFile does not exist');
                    end
                else
                    error ('MBDynInputFile is empty');
                end
            
                if isempty (options.CommMethod)
                    
                    error ('You must supply a communication method with ''CommMethod'' options if not supplying an mbsys object')
                else
                    
                    if ~isempty (options.SharedMemoryName)
                        self.commMethod = 'shared memory';
                        self.sharedMemoryName = options.SharedMemoryName;
                    elseif ~isempty (options.Path)
                        self.commMethod = 'local socket';
                        self.path = options.Path;
                    else
                        self.commMethod = 'inet socket';
                        self.host = options.Host;
                        self.port = options.Port;
                    end
                    
                end
                
                if (options.NNodes <= 0) && options.RefNode == 0
                
                    error ('MBCNodal:gaveportforlocal', ...
                        'NNodes and RefNode cannot both be zero');
                else
                    self.NNodes = options.NNodes;
                end

                % copy some of the chosen options over to the class properties
                self.useLabels = options.UseLabels;
                self.useAccelerations = options.UseAccelerations;
                self.useRefNode = options.UseRefNode;
                self.useMoments = options.UseMoments;
                self.nodeOrientiationType = options.NodeOrientationType;
            
            end
            
            if self.useLabels == true
                error ('MBCNodal:uselablesbug', ...
                    ['You have set UseLabels to true, unfortunately there is a bug in MBDyn that \n', ...
                     'means this option is not currently possible and will cause MBDyn to abort when \n', ...
                     'using external socket forces.']);
            end
                
            if isempty (options.OutputPrefix)
                [pathstr, name] = fileparts (self.MBDynInputFile);
                self.outputPrefix = fullfile (pathstr, name);
            else
                [pathstr, ~] = fileparts (options.OutputPrefix);
                if exist (pathstr, 'dir') ~= 7
                    error ('Output prefix directory does not exist');
                end
                self.outputPrefix = options.OutputPrefix;
            end
            
            if isempty (options.MBDynExecutable)
                self.MBDynExecutable = mbdyn.mint.find_mbdyn ();
            else
                if exist (options.MBDynExecutable, 'file') == 2
                    self.MBDynExecutable = options.MBDynExecutable;
                else
                    error ('The specified MBDynExecutable location does not exist');
                end
            end
            
            self.MBDynOutputFile = [self.outputPrefix, '.txt'];
            
            self.dataAndNext = options.DataAndNext;
            
        end

        function start (self, varargin)
            % Start mbdyn simulation
            %  
            % Syntax
            %
            % Start (mb)
            % Start (..., 'Parameter', Value)
            %
            % Input
            %
            % Some optional arguments may be supplied as parameter-value
            % pairs, where the required or desired options depend on the
            % particular configuration to be set up.
            %
            %
            % 'StartMBDyn' - determines whether the MBDyn process should be
            %   started. Default is true.
            %
            % 'Verbose' - 
            %
            % 'Timeout' - 
            %
            %
            
            options.Timeout = -1;
            options.Verbosity = 0;
            options.StartMBDyn = true;
            options.MBDynStartWaitTime = 1;
            
            options = parse_pv_pairs (options, varargin);
            
            mbdyn.pre.base.checkLogicalScalar (options.StartMBDyn, true, 'StartMBDyn');
            mbdyn.pre.base.checkNumericScalar (options.MBDynStartWaitTime, true, 'MBDynStartWaitTime');
            
            self.MBDynStartWaitTime = options.MBDynStartWaitTime;
            
            if options.StartMBDyn
                if options.Verbosity > 0
                    fprintf (1, 'Starting MBDyn\n');
                end
                % start mbdyn
                self.startMBdyn ('Verbosity', options.Verbosity);
            end
            
            % now initialise communication
            if strcmp (self.commMethod, 'local socket')
                
                try
                    self.cppcall ( 'Initialize', ...
                                   self.useRefNode, ...
                                   self.nodeOrientiationType, ...options.RefNodeOrientationType,
                                   self.NNodes, ...
                                   self.useLabels, ...
                                   self.nodeOrientiationType, ...
                                   self.useAccelerations, ...
                                   self.dataAndNext, ...
                                   options.Verbosity > 0, ...
                                   options.Timeout, ...
                                   'local', ...
                                   self.path );
                               
                catch err
                    % dump the mbdyn output to the command window (should
                    % be short at this stage
                    fprintf (1, 'MBDyn output:\n\n');
                    type (self.MBDynOutputFile);
                    rethrow (err);
                end
                           
                self.NNodes = GetNodes (self);
                
            elseif strcmp (self.commMethod, 'inet socket')
                
%                 if isempty (options.HostPort)
%                     error ('MBCNodal:noportforinet', ...
%                         'You have specified commethod ''inet'', but not specified a host port number with the ''HostPort'' option.');
%                 end
                try
                    self.cppcall ( 'Initialize', ...
                                   self.useRefNode, ...
                                   self.nodeOrientiationType, ...options.RefNodeOrientationType,
                                   self.NNodes, ...
                                   self.useLabels, ...
                                   self.nodeOrientiationType, ...
                                   self.useAccelerations, ...
                                   self.dataAndNext, ...
                                   options.Verbosity > 0, ...
                                   options.Timeout, ...
                                   'inet', ...
                                   self.host, ...
                                   self.port );
                               
                catch err
                    % dump the mbdyn output to the command window (should
                    % be short at this stage
                    fprintf (1, 'MBDyn output:\n\n');
                    type (self.MBDynOutputFile);
                    rethrow (err);
                end
                           
                self.NNodes = GetNodes (self);
                
            elseif strcmp (self.commMethod, 'shared memory')
                
%                 if isempty (options.HostPort)
%                     error ('MBCNodal:noportforinet', ...
%                         'You have specified commethod ''inet'', but not specified a host port number with the ''HostPort'' option.');
%                 end
                
                try
                    
                    self.cppcall ( 'Initialize', ...
                                   self.useRefNode, ...
                                   self.nodeOrientiationType, ...options.RefNodeOrientationType,
                                   self.NNodes, ...
                                   self.useLabels, ...
                                   self.nodeOrientiationType, ...
                                   self.useAccelerations, ...
                                   self.dataAndNext, ...
                                   options.Verbosity > 0, ...
                                   options.Timeout, ...
                                   self.sharedMemoryName );
                           
                catch err
                    % dump the mbdyn output to the command window (should
                    % be short at this stage
                    fprintf (1, 'MBDyn output:\n\n');
                    type (self.MBDynOutputFile);
                    rethrow (err);
                end
                           
                self.NNodes = GetNodes (self);
                
            else
                error ('MBCNodal:badcommmethod', ...
                    'Unrecognised communication method ''%s'' specified', commethod);
            end
        end
        
        function status = GetStatus (self)
            
            status = self.cppcall ( 'GetStatus');
            
        end

        function status = GetMotion (self)
            % obtain the last set of results from the mbdyn system
            %
            % Syntax
            %
            % GetMotion ()
            %
            % Description
            %
            % GetMotion makes available the last set of state data
            % generated by MBDyn. You must call GetMotion to get the first
            % set of data at the beginning of a simulation, and after each
            % call of 'applyForcesAndMoments' to get the latest data. Once
            % called the data can be accessed using the X, XP, XPP, and
            % similar methods.
            %
            % Input
            %
            %  None
            %
            % Output
            %
            %  status - zero if there are no errors, non-zero if there is a
            %    problem, or the simulation is complete, so no new data
            %    could be obtained.
            %
            
            status = self.cppcall ('GetMotion');
            
            self.needForces = true;
            self.needMoments = true;
            
        end

        function nnodes = GetNodes (self)
            % get the number of nodes in the system
            nnodes = self.cppcall ('GetNodes');
        end

        function label = KinematicsLabel (self, n)
            % gets the label associated with the 'n'th node
            
            if self.useLabels
                label = self.cppcall ('KinematicsLabel', n);
            else
                error ('MBCNodal:KinematicsLabel:nouselabels', ...
                    'You have set UseLabels to false, label data is not available.');
            end
        end

        function rot = GetRot (self)
            % get the rotation matrices for all nodes in the chosen format
            rot = self.cppcall ('GetRot');
            
            for ind = 1:self.NNodes
                
                if self.setNodePositions
                    
                    switch self.nodeOrientiationType
                        
                        case 'orientation matrix'
                            
                            om = mbdyn.pre.orientmat (self.nodeOrientiationType, rot(1:3,1:3,ind));
                            
                            
                        otherwise
                            
                            om = mbdyn.pre.orientmat (self.nodeOrientiationType, rot(1:3,ind));
                        
                    end
                    
                    self.structuralNodes{ind}.absoluteOrientation = om;
                    
                end
            
            end
        end
        
        function rot = GetRefNodeRot (self)
            % get the rotation matrix in the chosen format
            rot = self.cppcall ('GetRefNodeRot');
        end
        
        function pos = NodePositions (self, n)
            % gets the positions of one or more nodes
            %
            % Syntax
            %
            % NodePositions ()
            % NodePositions (n)
            %
            % Input
            %
            %  n - vector of one or more node numbers for which to get the
            %    position. If not supplied, the positions of all nodes will
            %    be returned.
            %
            % Output
            %
            %  pos - (3 x k) matrix of k node positions, one for each
            %    node number supplied in input 'n'. Each column represents
            %    a node. Will contain the positions of all nodes if n is
            %    not supplied
            %
            %
            
            if nargin < 2
                n = 1:self.NNodes;
            end
            
            pos = nan * ones (3, numel(n));
            
            for ind = 1:numel (n)
                
               pos(1:3,ind) = X (self, n(ind));
               
               if self.setNodePositions
                    self.structuralNodes{n(ind)}.absolutePosition = pos(1:3,ind);
               end
               
            end
            
        end
        
        function vel = NodeVelocities (self, n)
            % gets the velocities of one or more nodes
            %
            % Syntax
            %
            % NodeVelocities ()
            % NodeVelocities (n)
            %
            % Input
            %
            %  n - vector of one or more node numbers for which to get the
            %    velocity. If not supplied, the velocities of all nodes
            %    will be returned.
            %
            % Output
            %
            %  vel - (3 x k) matrix of k node velocities, one for each
            %    node number supplied in input 'n'. Each column represents
            %    a node. Will contain the velocities of all nodes if n is
            %    not supplied
            %
            %
            
            if nargin < 2
                n = 1:self.NNodes;
            end
            
            vel = nan * ones (3, numel(n));
            
            for ind = 1:numel (n)
                
               vel(1:3,ind) = XP (self, n(ind));
               
               if self.setNodePositions
                    self.structuralNodes{n(ind)}.absoluteVelocity = vel (1:3,ind);
               end
               
            end
            
        end
        
        function accel = NodeAccelerations (self, n)
            % gets the accelerations of one or more nodes
            %
            % Syntax
            %
            % NodeAccelerations ()
            % NodeAccelerations (n)
            %
            % Input
            %
            %  n - vector of one or more node numbers for which to get the
            %    acceleration. If not supplied, the velocities of all nodes
            %    will be returned.
            %
            % Output
            %
            %  accel - (3 x k) matrix of k node accelerations, one for each
            %    node number supplied in input 'n'. Each column represents
            %    a node. Will contain the accelerations of all nodes if n
            %    is not supplied
            %
            %
            
            if nargin < 2
                n = 1:self.NNodes;
            end
            
            accel = nan * ones (3, numel(n));
            
            for ind = 1:numel (n)
                
               accel(1:3,ind) = XPP (self, n(ind));
               
            end
            
        end
        
        function theta = NodeThetas (self, n)
            % gets the angular positions of one or more nodes
            %
            % Syntax
            %
            % NodeThetas ()
            % NodeThetas (n)
            %
            % Input
            %
            %  n - vector of one or more node numbers for which to get the
            %    angular position. If not supplied, the angular positions
            %    of all nodes will be returned.
            %
            % Output
            %
            %  pos - (3 x k) matrix of k node angular positions, one for each
            %    node number supplied in input 'n'. Each column represents
            %    a node. Will contain the angular positions of all nodes if
            %    n is not supplied
            %
            %
            
            if nargin < 2
                n = 1:self.NNodes;
            end
            
            theta = nan * ones (3, numel(n));
            
            for ind = 1:numel (n)
                
               theta(1:3,ind) = Theta (self, n(ind));
               
               if self.setNodePositions
                   
                    om = mbdyn.pre.orientmat (self.nodeOrientiationType, theta(1:3,ind));
                   
                    self.structuralNodes{n(ind)}.absoluteOrientation = om.orientationMatrix;
               end
               
            end
            
        end
        
        function vel = NodeOmegas (self, n)
            % gets the angular velocities of one or more nodes
            %
            % Syntax
            %
            % NodeOmegas ()
            % NodeOmegas (n)
            %
            % Input
            %
            %  n - vector of one or more node numbers for which to get the
            %    angular velocity. If not supplied, the angular velocities
            %    of all nodes will be returned.
            %
            % Output
            %
            %  vel - (3 x k) matrix of k node angular velocities, one for
            %    each node number supplied in input 'n'. Each column
            %    represents a node. Will contain the angular velocities of
            %    all nodes if n is not supplied
            %
            %
            
            if nargin < 2
                n = 1:self.NNodes;
            end
            
            vel = nan * ones (3, numel(n));
            
            for ind = 1:numel (n)
                
               vel (1:3,ind) = Omega (self, n(ind));
               
               if self.setNodePositions
                    self.structuralNodes{n(ind)}.absoluteAngularVelocity = vel (1:3,ind);
               end
               
            end
            
        end
        
        function accel = NodeAngularAccels (self, n)
            % gets the angular accelerations of one or more nodes
            %
            % Syntax
            %
            % NodeAngularAccels ()
            % NodeAngularAccels (n)
            %
            % Input
            %
            %  n - vector of one or more node numbers for which to get the
            %    angular acceleration. If not supplied, the angular
            %    accelerations of all nodes will be returned.
            %
            % Output
            %
            %  accel - (3 x k) matrix of k node angular accelerations, one
            %    for each node number supplied in input 'n'. Each column
            %    represents a node. Will contain the angular accelerations
            %    of all nodes if n is not supplied
            %
            %
            
            if nargin < 2
                n = 1:self.NNodes;
            end
            
            accel = nan * ones (3, numel(n));
            
            for ind = 1:numel (n)
                
               accel (1:3,ind) = OmegaP (self, n(ind));
               
            end
            
        end
        
        function pos = X (self, n)
            % gets the position of a single node with number n
            %
            % Syntax
            %
            %  X (n)
            %
            % Input
            %
            %  n - scalar integer representing the node number for which
            %    the position is to be returned
            %
            % Output
            %
            %  pos - (3 x 1) vector containing the xyz position of the node
            %
            pos = self.cppcall ('X', n)';
        end

        function vel = XP (self, n)
            % gets the velocity of a single node with number n
            %
            % Syntax
            %
            %  XP (n)
            %
            % Input
            %
            %  n - scalar integer representing the node number for which
            %    the velocity is to be returned
            %
            % Output
            %
            %  pos - (3 x 1) vector containing the xyz velocity of the node
            %
            vel = self.cppcall ('XP', n)';
        end
        
        function acc = XPP (self, n)
            % gets the acceleration of a single node with number n
            %
            % Syntax
            %
            %  XPP (n)
            %
            % Input
            %
            %  n - scalar integer representing the node number for which
            %    the acceleration is to be returned
            %
            % Output
            %
            %  pos - (3 x 1) vector containing the xyz acceleration of the
            %    node
            %
            
            if ~self.useAccelerations
                error ('MBCNodal:xpp:nouseaccelerations', ...
                    'You have set UseAccelerations to false, acceleration data is not available.')
            else
                acc = self.cppcall ('XPP', n)';
            end
        end
        
        function theta = Theta (self, n)
            % gets the angular position of a single node with number n
            %
            % Syntax
            %
            %  Theta (n)
            %
            % Input
            %
            %  n - scalar integer representing the node number for which
            %    the angular position is to be returned
            %
            % Output
            %
            %  pos - (3 x 1) vector containing the xyz angular position of
            %    the node
            %
            theta = self.cppcall ('Theta', n)';
        end
        
        function theta = Euler123 (self, n)
            % gets the angular position of a single node with number n
            %
            % Syntax
            %
            %  Theta (n)
            %
            % Input
            %
            %  n - scalar integer representing the node number for which
            %    the angular position is to be returned
            %
            % Output
            %
            %  pos - (3 x 1) vector containing the xyz angular position of
            %    the node
            %
            theta = self.cppcall ('Euler123', n)';
        end
        
        function w = Omega (self, n)
            % gets the angular velocity of a single node with number n
            %
            % Syntax
            %
            %  Omega (n)
            %
            % Input
            %
            %  n - scalar integer representing the node number for which
            %    the angular velocity is to be returned
            %
            % Output
            %
            %  pos - (3 x 1) vector containing the xyz angular velocity of
            %    the node
            %
            w = self.cppcall ('Omega', n)';
        end
        
        function w = OmegaP (self, n)
            % gets the angular acceleration of a single node with number n
            %
            % Syntax
            %
            %  OmegaP (n)
            %
            % Input
            %
            %  n - scalar integer representing the node number for which
            %    the angular acceleration is to be returned
            %
            % Output
            %
            %  pos - (3 x 1) vector containing the xyz anangular
            %    acceleration of the node
            %
            if ~self.useAccelerations
                error ('MBCNodal:omegap:nouseaccelerations', ...
                    'You have set UseAccelerations to false, angular acceleration data is not available.');
            else
                w = self.cppcall ('OmegaP', n)';
            end
        end
        
        function F (self, forces)
            % sets the nodal forces (for all structural external nodes)
            %
            % Syntax
            %
            % F (forces)
            %
            % Input
            %
            %  forces - (3 x n) matrix of forces in the x, y and z
            %    directions for all n nodes, one node for each column
            %
            % Output
            %
            %  None
            %
            
            self.cppcall ('F', forces);
            self.needForces = false;
            
            if ~self.useMoments
                self.cppcall ('M', zeros (size (forces)));
                self.needMoments = false;
            end
        end
        
        function M (self, moments)
            % sets the nodal moments (for all structural external nodes)
            %
            % Syntax
            %
            % M (forces)
            %
            % Input
            %
            %  moments - (3 x n) matrix of forces about the x, y and z
            %    axes for all n nodes, one node for each column
            %
            % Output
            %
            %  None
            %
            
            if self.useMoments
                self.cppcall ('M', moments );
                self.needMoments = false;
            else
                warning ('MBCNodal:nousemoments', ...
                    'You have set UseMoments to false, but attempted to apply moments, these will be ignored.')
            end
        end
        
        function result = applyForcesAndMoments (self, convergence_flag)
            % sends the forces and moments to the mbdyn system 
            %
            % Syntax
            %
            % result = applyForcesAndMoments (self, convergence_flag)
            %
            % Input
            %
            %  convergence_flag - boolean flag indicating whether the
            %    external system has converged. When this flag is false,
            %    MBDyn will no advance the time step, but rather
            %    recalculate the state of the system and return a new set
            %    of state variables (positions, velocities etc.) based on
            %    the new forces. When the flag is true, it indicates the
            %    system has converged, this is the final set of forces for
            %    this time step, and MBDyn will apply the forces and
            %    advance to the next time step using these forces.
            %
            % Output
            %
            %  result - 
            
            if self.needForces
                error ('MBCNodal:notsetforces', ...
                    'You must set the nodal forces in the system before applying them');
            end
            
            if self.needMoments
                error ('MBCNodal:notsetforces', ...
                    'You must set the nodal moments in the system before applying them');
            end
            
            % Sends the forces to the mbdyn system
            result = self.cppcall ('PutForces', logical (convergence_flag) );
        end

    end
    
    methods (Access = protected)
        
        function startMBdyn (self, varargin)
            % run mbdyn with the appropriate commands
            
            options.Verbosity = 0;
            
            options = parse_pv_pairs (options, varargin);
            
            if options.Verbosity > 0
                Pcmds = ['-', repmat('P', [1, options.Verbosity])];
            else
                Pcmds = '';
            end
            
            % start mbdyn
            cmdline = sprintf ('%s %s -f "%s" -o "%s" > "%s" 2>&1 &', ...
                                self.MBDynExecutable, ...
                                Pcmds, ...
                                self.MBDynInputFile, ...
                                self.outputPrefix, ...
                                self.MBDynOutputFile  ...
                                                 );
                                             
            if options.Verbosity > 0
                fprintf (1, 'Starting MBDyn with command:\n%s\n', cmdline);
            end

            [status, cmdout] = cleansystem ( cmdline );
            
            pause (self.MBDynStartWaitTime);
            
        end
        
        
        
        
        function status = cppstatus2statusenum (cppstatus)
        
            switch cppstatus
                
                case -1
                    
                case -2
                    
                case 0
                    
                case 1
                    
                otherwise
                    
                    error ('Unknown status number');
                    
            end
            
        end
             
        
    end

end
