classdef MBCNodal < cppinterface
    % MBCNodal: interface for MBDyn multibody dynamics library using sockets. 
    %
    % Description
    %
    % 

    properties (SetAccess = private, GetAccess = public)
        
        NNodes; % Number of structural external force nodes
        
        useLabels; % Flag indicating whether to make node label numbers available
        useAccelerations; % Flag indicating whether accelerations are to be made available
        useRefNode; % Flag indicating whether there is a reference node
        
        useMoments; % Flag indicating whether moments are to be used
        
    end
    
    properties (SetAccess = private, GetAccess = private)
        
        % needMoments -Internal flag used to determine if we need to apply
        % moments before advancing
        needMoments;
        
        % needForces -Internal flag used to determine if we need to apply
        % forces before advancing 
        needForces;
        
    end
    
    methods
        
        % Constructor
        function self = MBCNodal ()

            % initialise the cppinterface parent class by passing the
            % mexfunction to the superclass constructor
            self = self@cppinterface(@mbdyn.mexMBCNodal);

        end

        function Initialize (self, commethod, comstring, varargin)
            %  Initialize nodal mbdyn simulation
            %  
            % Syntax
            %
            % MBCNodal.Initialize (commethod, commpath)
            % MBCNodal.Initialize (..., 'Parameter', Value)
            %
            % Input
            %
            %  commethod - string determining what type of socket
            %    communication to use to communicate with MBDyn. The
            %    avaialble options are: 'local' and 'inet'. The 'local'
            %    options used a file based method and you must supply the
            %    file location in the commpath argument (see below). The
            %    'inet' method used TCP/IP over a network and the address
            %    to use must be supplied in commpath. If using 'inet' you
            %    must also supply a host port to use for communication
            %    using the 'HostPort' optional argument described below.
            %
            %  comstring - a string containing the communication path to
            %    use. The contents depends on the value of commethod above.
            %
            % Further arguments may be supplied as parameter-value pairs,
            % where the required or desired options depend on the
            % particular configuration to be set up.
            %
            % 'HostPort' - if the inet communication method is to be used,
            %   this is the port number to use for communication, and must
            %   be supplied.
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
            
            options.UseRefNode = false; 
            options.RefNodeOrientationType = 'none'; % refnode rotation type 'none'. Will be ignored if UseRefNode is false
            options.NNodes = 0; % number of nodes, should match up with problem file?
            options.NodeOrientationType = 'mat';
            options.UseLabels = true;
            options.UseAccelerations = false; % don't handle accelerations by default
            options.HostPort = []; % inet port to use with commethod 'inet'
            options.UseMoments = false;
            options.DataAndNext = true;
            options.Timeout = -1;
            options.Verbose = false;
            
            options = parse_pv_pairs (options, varargin);
            
            if (options.NNodes <= 0) && options.RefNode == 0
                
                error ('MBCNodal:gaveportforlocal', ...
                    'NNodes and RefNode cannot both be zero');
                
            end
            
            if options.UseLabels == true
                error ('MBCNodal:uselablesbug', ...
                    ['You have set UseLabels to true, unfortunately there is a bug in MBDyn that \n', ...
                     'means this option is not currently possible and will cause MBDyn to abort when \n', ...
                     'using external socket forces.']);
            end
            
            % copy some of the chosen options over to the class properties
            self.useLabels = options.UseLabels;
            self.useAccelerations = options.UseAccelerations;
            self.useRefNode = options.UseRefNode;  
            
            self.useMoments = options.UseMoments;
            
            if strcmp (commethod, 'local')
                
                if ~isempty (options.HostPort)
                    warning ('MBCNodal:gaveportforlocal', ...
                        'You have specified a port number with the ''local'' connection method which is file based.');
                end
                
                self.cppcall ( 'Initialize', ...
                               options.UseRefNode, ...
                               options.RefNodeOrientationType, ...
                               options.NNodes, ...
                               options.UseLabels, ...
                               options.NodeOrientationType, ...
                               options.UseAccelerations, ...
                               options.DataAndNext, ...
                               options.Verbose, ...
                               options.Timeout, ...
                               commethod, ...
                               comstring );
                           
                self.NNodes = GetNodes (self);
                
            elseif strcmp (commethod, 'inet')
                
                if isempty (options.HostPort)
                    error ('MBCNodal:noportforinet', ...
                        'You have specified commethod ''inet'', but not specified a host port number with the ''HostPort'' option.');
                end
                
                self.cppcall ( 'Initialize', ...
                               options.UseRefNode, ...
                               options.RefNodeOrientationType, ...
                               options.NNodes, ...
                               options.UseLabels, ...
                               options.NodeOrientationType, ...
                               options.UseAccelerations, ...
                               options.DataAndNext, ...
                               options.Verbose, ...
                               options.Timeout, ...
                               commethod, ...
                               comstring, ...
                               options.HostPort );
                           
                self.NNodes = GetNodes (self);
                
            else
                error ('MBCNodal:badcommmethod', ...
                    'Unrecognised communication method ''%s'' specified', commethod);
            end
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
                
               pos (1:3,ind) = X(self, n(ind));
               
            end
            
        end
        
        function pos = X(self, n)
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
            result = self.cppcall ('PutForces', boolean (convergence_flag) );
        end

    end

end
