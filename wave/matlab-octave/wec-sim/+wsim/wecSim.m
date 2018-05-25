classdef wecSim < handle
    
    properties
        wavePlotNumPointsX = 50;
        wavePlotNumPointsY = 50;
        wavePlotDomainSize = [];
    end
    
    properties (GetAccess = public, SetAccess = private)
        
        powerTakeOffs;
        forceRelTol;
        forceAbsTol;
        mBDynOutputFile;
        mBDynInputFile;
        caseDirectory;

        loggingSettings;
        
    end
    
    properties (GetAccess = private, SetAccess = private)
        
        ptoIndexMap;
        hydroNodeIndexMap;
        mBDynSystem;
        hydroSystem;
        readyToRun;
        simComplete;
        hydroMotionSyncStep = 1;
        hydroMotionSyncStepCount = 1;
        
        
        logger;
        nMBDynNodes;
        mBDynPostProc;
        
        % logging variables
        lastTime;                   % last value of simulation time
        lastPositions;              % last calculated cartesian positions
        lastAngularPositions;       % last calculated angular positions (euler angles)
        lastVelocities;             % last calculated velocites
        lastAngularVelocities;      % last calculated angular velocities
        lastAccelerations;          % last calculated accelerations
        lastAngularAccelerations;   % last calculated angular accelerations
        lastForceHydro;             % last calculated total of all hydrodynamic forces
        lastForceExcitation;        % last calculated value of all hydrodynamic excitation forces
        lastForceExcitationRamp;    % last calculated value of all hydrodynamic excitation forces with applied ramp
        lastForceExcitationLin;     % last calculated value of hydrodynamic linear excitation forces
        lastForceExcitationNonLin;  % last calculated value of hydrodynamic nonlinear excitation forces
        lastForceRadiationDamping;  % last calculated value of hydrodynamic radiation damping forces
        lastForceRestoring;         % last calculated value of hydrodynamic restoring (buoyancy) forces
        lastForceMorrison;          % last calculated value of hydrodynamic morrison forces
        lastForceViscousDamping;    % last calculated value of hydrodynamic viscous damping forces
        lastForceAddedMassUncorrected; % last calculated value of hydrodynamic forces due to added mass (uncorrected)
        lastNodeForcesAndMomentsUncorrected;
        
        wavePlotX;
        wavePlotY;
        
        drawAxesH;
        waveSurfH;
        
    end
    
    methods
        
        function self = wecSim (hsys, mbsys, varargin)
            % wsim.wecSim constructor
            %
            
            options.PTOs = {};
%             options.StartTime = [];
%             options.EndTime = [];
            options.NEMOHSim = [];
            options.LoggingSettings = wsim.loggingSettings ();
            options.MBDynVerbosity = 0;
            options.OutputDir = '';
            
            options = parse_pv_pairs (options, varargin);
            
            assert (isa (hsys, 'wsim.hydroSystem'), ...
                'hsys must be an wsim.hydroSystem object');
            
            assert (isa (mbsys, 'mbdyn.pre.system'), ...
                'mbsys must be an mbdyn.pre.system object');
            
            if ~isempty (options.NEMOHSim)
                assert (isa (options.NEMOHSim, 'nemoh.simulation'), ...
                    'NEMOHSim must be a nemoh.simulation object' ); 
            end
            
            self.caseDirectory = hsys.simu.caseDir;
            self.loggingSettings = options.LoggingSettings;
            self.mBDynSystem = mbsys;
            self.hydroSystem = hsys;
            self.readyToRun = false;
            self.simComplete = false;
            self.nMBDynNodes = nan;
            
            % add the spplies PTOs
            if ~isempty (options.PTOs)
                if ~iscell (options.PTOs)
                    options.PTOs = {options.PTOs};
                end
                for ptoind = 1:numel (options.PTOs)
                    self.addPTO (options.PTOs{ptoind});
                end
            end

        end
        
        function addPTO (self, pto)
            
            assert (isa (pto, 'wsim.powerTakeOff'), ...
                'pto must be a wsim.powerTakeOff object (or derived class)');
            
            % check that the PTO nodes are in the MBDyn system external
            % structual force element
            strinfo = self.mBDynSystem.externalStructuralInfo ();
            
            for ptoind = 1:numel(self.powerTakeOffs)
                
                foundptonode = false;
                for nodeind = 1:numel (strinfo.Nodes)
                    
                    if strinfo.Nodes{nodeind} == self.powerTakeOffs{ptoind}.referenceNode
                        
                        foundptonode = true;
                        
                        break;
                    
                    end
                
                end
                
                if foundptonode == false
                    error ('Could not find reference node for new PTO in MBDyn system structural external force element');
                end
                
                foundptonode = false;
                for nodeind = 1:numel (strinfo.Nodes)
                    
                    if strinfo.Nodes{nodeind} == self.powerTakeOffs{ptoind}.otherNode
                        
                        foundptonode = true;
                    
                        break;
                        
                    end
                
                end
                
                if foundptonode == false
                    error ('Could not find non-reference node for new PTO in MBDyn system structural external force element');
                end
                
            end
            
            % append it to the existing PTO objects
            self.powerTakeOffs = [ self.powerTakeOffs, {pto}];
            
            % set the id of the pto
            self.powerTakeOffs{end}.id = numel (self.powerTakeOffs);
            
            % mark ready to run false, as the PTO index map will need to be
            % updated before proceeding to run a simulation
            self.readyToRun = false;
            
        end
        
        function prepare (self)
            
            self.mapPTOForceInds ();
            self.mapHydroForceInds ();
            self.initDataStructures ();
            
            % TODO: check if NEMOH data needs updated
            
            % ensure hydro system transient simulation is ready
%             self.hydroSystem.timeDomainSimSetup ();

            % clear any mbdyn.postproc object
            self.mBDynPostProc = [];
            
            self.readyToRun = true;
            self.simComplete = false;
            
        end
        
        function [ datalog, mbdyn_postproc ] = run (self, varargin)
            % Run the WEC simulation
            
            options.MBDynInputFile = '';
            options.OutputFilePrefix = fullfile (self.caseDirectory, ['output_', datestr(now (), 'yyyy-mm-dd_HH-MM-SS-FFF')]);
            options.Verbosity = 0;
            options.AbsForceTolerance = 100;
            options.RelForceTolerance = 1e-5;
            options.MinIterations = 0;
            options.MaxIterations = self.mBDynSystem.problems{1}.maxIterations;
            options.TimeExecution = false;
            options.HydroMotionSyncSteps = 1;
            options.ForceMBDynNetCDF = true;
            
            if isempty (options.MBDynInputFile)
                self.mBDynInputFile = fullfile (options.OutputFilePrefix, 'mbdyn_input_file.mbd');
            else
                assert (ischar (options.MBDynInputFile), ...
                    'MBDynInputFile must be string containing the file path where the MBDyn input file should be generated');
            end
            
            options = parse_pv_pairs (options, varargin);
            
            check.isPositiveScalarInteger (options.Verbosity, true, 'Verbosity', true);
            
            assert (ischar (options.OutputFilePrefix), ...
                'OutputFilePrefix must be a string');
            
            if exist (options.OutputFilePrefix, 'dir') == 0
                mkdir (options.OutputFilePrefix);
            end
                
            % check for positive numeric scalar
            check.isNumericScalar (options.AbsForceTolerance, true, 'AbsForceTolerance', 1);
            check.isNumericScalar (options.RelForceTolerance, true, 'RelForceTolerance', 1);
            
            check.isPositiveScalarInteger (options.MaxIterations, true, 'Verbosity');
            
            assert (options.MinIterations <= options.MaxIterations, ...
                'MinIterations (%d) is not less than or equal to MaxIterations (%d)', ...
                options.MinIterations,  options.MaxIterations)
            
            check.isLogicalScalar (options.TimeExecution, true, 'TimeExecution');
            
            check.isPositiveScalarInteger (options.HydroMotionSyncSteps, true, 'HydroMotionSyncSteps');
            
            check.isLogicalScalar (options.ForceMBDynNetCDF, true, 'ForceMBDynNetCDF');
            
            % -----------------   input checking finished
            
            if self.readyToRun == false
                error ('Simulation is not ready to run, have you run ''prepare'' yet?');
            end
            
            if nargout > 1 || options.ForceMBDynNetCDF
                % ensure mbdyn will output a netcdf file so we can load the
                % results from this after the sim
                if isempty (self.mBDynSystem.controlData.OutputResults)
                    self.mBDynSystem.controlData.OutputResults = {'netcdf', 'no text'};
                end
            end
            
            self.hydroMotionSyncStep = options.HydroMotionSyncSteps;
            % set the hydro--mbdyn motion sync step count equal to the
            % number of steps to sync on so that forces are calculated the
            % first time applyHydroForces is called (where this step count
            % is used)
            self.hydroMotionSyncStepCount = self.hydroMotionSyncStep;
            
            siminfo.TStart = self.mBDynSystem.problems{1}.initialTime;
            siminfo.TEnd = self.mBDynSystem.problems{1}.finalTime;
            siminfo.TStep = self.mBDynSystem.problems{1}.timeStep;
            siminfo.MBDynSystem = self.mBDynSystem;
            siminfo.HydroSystem = self.hydroSystem;
            siminfo.HydroMotionSyncSteps = options.HydroMotionSyncSteps;
            siminfo.OutputDirectory = options.OutputFilePrefix;
            siminfo.CaseDirectory = self.caseDirectory;
            
            % start the PTO components
            self.startPTOs (siminfo);
            
            % generate input file and start mbdyn

            % create the communicator object. As an mbdyn system object is
            % supplied, the mbdyn input file will be generated
            % automatically
            mb = mbdyn.mint.MBCNodal ('MBDynPreProc', self.mBDynSystem, ...
                'UseMoments', true, ...
                'MBDynInputFile', self.mBDynInputFile, ...
                'OverwriteInputFile', true, ...
                'OutputPrefix', options.OutputFilePrefix, ...
                'NodeOrientationType', 'euler 123' ...
                );
            
            % copy over the input file location to make it easier to
            % examine later if required
            self.mBDynOutputFile = mb.MBDynOutputFile;
            
%             % ensure MBCNodal is destroyed in the event of a problem
%             % (closes communication to MBDyn and tells it to quit so
%             % sockets and so on are also cleaned up)
%             if ~isoctave
%                 % FIXME: put onCleanup back when Octave is compatible
%                 CC = onCleanup (@() delete (mb));
%             end
            
            mb.start ('Verbosity', options.Verbosity);
            
            % get the number of nodes in the problem
            self.nMBDynNodes = mb.GetNodes ();
            
            % preallocate variables to hold total forces and moments at
            % each time step
            forces_and_moments = zeros (6, self.nMBDynNodes);
            
            % take initial time from the MBDyn system problem
            self.lastTime = siminfo.TStart;
            
            % fetch the initial configuration of the nodes from MBDyn and
            % check for errors
            status = mb.GetMotion ();
            
            if status ~= 0
                self.readyToRun = false;
                if exist (self.mBDynOutputFile, 'file')
                    self.displayLastNLinesOfFile (self.mBDynOutputFile, 50);
                end
                error ('mbdyn returned %d, aborting sim, check output file:\n%s\nfor clues at to why this happened.', ...
                        status, mb.MBDynOutputFile)
            end
            
            % get the current motion of the multibody system which was sent
            % by MBDyn
            [pos, vel, accel] = self.getMotion (mb);
                  
            % calculate and apply hydrodynamic forces based on the motion
            forces_and_moments = self.applyHydroForces (forces_and_moments, self.lastTime, pos, vel, accel);
            
            % get and apply the PTO forces and moments, in this case the
            % position and velocities etc are obtained from the mbsys
            % object stored in the mb object, hence no need to supply them.
            % The node positions etc. are updated by the call to
            % mb.GetMotion (). This is possible because all PTO objects
            % also have access to the same mbsys object
            forces_and_moments = self.applyPTOForces (forces_and_moments);

            % set the new forces an moments ready to be sent to MBDyn
            mb.F (forces_and_moments(1:3,:));
            mb.M (forces_and_moments(4:6,:));
            
            % send the forces and moments to MBDyn, but noting that we have
            % not yet converged (so MBDyn will send new motion based on
            % these forces, and not advance the self.lastTime step)
            mbconv = mb.applyForcesAndMoments (false);
            
            % store most recently calculated motion so it can be logged
            self.lastPositions = pos(1:3,:);
            self.lastAngularPositions = pos(4:6,:);
            self.lastVelocities = vel(1:3,:);
            self.lastAngularVelocities = vel(4:6,:);
            self.lastAccelerations = accel(1:3,:);
            self.lastAngularAccelerations = accel(4:6,:);
            self.lastNodeForcesAndMomentsUncorrected = forces_and_moments;

            self.advanceStep ();
            
            % now begin the simulation loop, beginning from the next time
            % index
            ind = 2;
            if options.TimeExecution, tic; end
            while status == 0
                
                status = mb.GetMotion ();
                
                if status ~= 0
                    self.readyToRun = false;
                    continue;
                end
                
                self.lastTime = self.lastTime + self.mBDynSystem.problems{1}.timeStep;

                % get the current motion of the multibody system
                [pos, vel, accel] = self.getMotion (mb);
                    
                % clear out the previous forces and moments
                forces_and_moments = zeros (6, self.nMBDynNodes);

                % calculate new hydrodynamic interaction forces (also
                % updates self.lastForceHydro with new values)
                forces_and_moments = self.applyHydroForces (forces_and_moments, self.lastTime, pos, vel, accel);

                % PTO forces
                forces_and_moments = self.applyPTOForces (forces_and_moments);

                mb.F (forces_and_moments(1:3,:));
                mb.M (forces_and_moments(4:6,:));

                mbconv = mb.applyForcesAndMoments (false);

                status = mb.GetMotion ();

                if status ~= 0
                    % quit, MBDyn has finished simulation (or some error
                    % occured)
                    self.readyToRun = false;
                    ind = ind + 1;

                    continue;
                end

                % repeat the force calulation to test convergence

                % clear out the previous forces and moments
                forces_and_moments = zeros (6, self.nMBDynNodes);

                % get the current motion of the multibody system
                [pos, vel, accel] = self.getMotion (mb);

                % get the previous calculated values of hydrodynamic forces
                hydroforces = self.lastForceHydro;

                % calculate new hydrodynamic interaction forces (also
                % updates self.lastForceHydro with new values)
                forces_and_moments = self.applyHydroForces (forces_and_moments, self.lastTime, pos, vel, accel);

                % PTO forces
                forces_and_moments = self.applyPTOForces (forces_and_moments);

                mb.F (forces_and_moments(1:3,:));
                mb.M (forces_and_moments(4:6,:));

                mbconv = mb.applyForcesAndMoments (false);

                % check for force convergence
                forcediff = abs (hydroforces - self.lastForceHydro);

                maxforces = max(hydroforces, self.lastForceHydro);
                relforcediff = abs(forcediff) ./ abs(maxforces);
                relforcediff(maxforces == 0) = 0;
                itercount = 1;

                % iterate until force/motion convergence (or max
                % iterations)
                while mbconv ~= 0 ...
                        || itercount < options.MinIterations ...
                        || (max (forcediff(:)) > options.AbsForceTolerance) ...
                        || (ind > 3 && (max (relforcediff(:)) > options.RelForceTolerance))

                    status = mb.GetMotion ();

                    if status ~= 0
                        self.readyToRun = false;
                        break;
                    end

                    % clear out the previous forces and moments
                    forces_and_moments = zeros (6, self.nMBDynNodes);

                    % get the current motion of the multibody system
                    [pos, vel, accel] = self.getMotion (mb);

                    % get the previous calculated values of hydrodynamic forces
                    hydroforces = self.lastForceHydro;

                    % calculate new hydrodynamic interaction forces (also
                    % updates self.lastForceHydro with new values)
                    forces_and_moments = self.applyHydroForces (forces_and_moments, self.lastTime, pos, vel, accel);

                    % PTO forces
                    forces_and_moments = self.applyPTOForces (forces_and_moments);

                    mb.F (forces_and_moments(1:3,:));
                    mb.M (forces_and_moments(4:6,:));

                    mbconv = mb.applyForcesAndMoments (false);

                    itercount = itercount + 1;

                    if itercount > options.MaxIterations
                        % TODO: exit gracefully if max iterations exceeded
                        error ('mbdyn iterations exceeded max allowed');
                    end

                    forcediff = abs (hydroforces - self.lastForceHydro);
                    maxforces = max(hydroforces, self.lastForceHydro);
                    relforcediff = abs(forcediff) ./ abs(maxforces);
                    relforcediff(maxforces == 0) = 0;

                end


                % get latest motion from MBDyn
                status = mb.GetMotion ();
                
                if status ~= 0
                    
                    % quit, MBDyn sim is finished or there's an error
                    self.readyToRun = false;
                    
                    ind = ind + 1;
                    
                    break;
                end
                
                % apply the last calculated forces from the iteration to
                % the system
                mb.F (forces_and_moments(1:3,:));
                mb.M (forces_and_moments(4:6,:));
                
                mbconv = mb.applyForcesAndMoments (true);
                
                % store most recently calculated motion so it can be logged
                self.lastPositions = pos(1:3,:);
                self.lastAngularPositions = pos(4:6,:);
                self.lastVelocities = vel(1:3,:);
                self.lastAngularVelocities = vel(4:6,:);
                self.lastAccelerations = accel(1:3,:);
                self.lastAngularAccelerations = accel(4:6,:);
                self.lastNodeForcesAndMomentsUncorrected = forces_and_moments;
                
                self.advanceStep ();
                
                ind = ind + 1;
                
            end
            
            fprintf ( 1, 'Reached time %fs (%fs of an intended %fs), in %d steps, postprocessing ...\n', ...
                      self.lastTime, ...
                      self.lastTime - siminfo.TStart, ...
                      siminfo.TEnd - siminfo.TStart, ...
                      ind-1 );
            
            if self.loggingSettings.nodeForcesAndMoments || self.loggingSettings.forceAddedMass
                
                self.logger.setSeries ('NodeForcesAndMoments', self.logger.data.NodeForcesAndMomentsUncorrected);

                % the total forces on hydrodynamic bodies must be corrected
                [ self.logger.data.NodeForcesAndMoments(:,self.hydroNodeIndexMap(:,1),:), ...
                  F_added_mass ] = ...
                    correctAddedMassForce ( self.hydroSystem, ...
                                            self.logger.data.NodeForcesAndMomentsUncorrected(:,self.hydroNodeIndexMap(:,1),:), ...
                                            self.logger.data.ForceAddedMassUncorrected, ...
                                            [ self.logger.data.Accelerations; self.logger.data.AngularAccelerations ] );

                if self.loggingSettings.forceAddedMass
                    self.logger.setSeries ('ForceAddedMass', F_added_mass);
                end
            
            end
            
            % tell the PTOs we are done
            self.finishPTOs ();
            
            if options.TimeExecution, toc; end
            
            fprintf (1, 'Simulation complete\n');
            
            self.readyToRun = false;
            
            self.logger.truncateAllVariables ();
            
            datalog = self.logger;
            
            self.simComplete = true;
            
            if nargout > 1
                mbdyn_postproc = mbdyn.postproc (self.mBDynInputFile, self.mBDynSystem);
                self.mBDynPostProc = mbdyn_postproc;
            end
            
        end
        
        
        function plotdata = drawStep (self, tind, varargin)
            % draw the system at the given time step index
            %
            % Syntax
            %
            % plotdata = drawStep (wsobj, tind, 'Parameter', value)
            %
            % Input
            %
            %  wsobj - wsim.wecSim object
            %
            %  Addtional optional arguments are supplied as parameter-value
            %  pairs. The available options are:
            %
            %  'PlotAxes' - handle to plot axes in which to draw the
            %    trajectories. If not supplied a new figure and axes are
            %    created for the plot.
            %
            %  'AxLims' - (3 x 2) matrix containing limits the plot axes
            %    should be set to, each row is the x, y and z axes limits
            %    respectively. Uselful to focus on a particular region of
            %    the system. If not supplied, suitable axes limits will be
            %    determined based on the motion data.
            %
            %  'Title' - flag determining whether to add a title to the
            %    system plot. Default is true.
            %
            %  'DrawMode' - string indicating the drawing mode. Can be one
            %    of: 'solid', 'ghost', 'wire', 'wireghost'.
            %
            %  'DrawLabels' - flag indicating whether to draw node labels,
            %    default is false.
            %
            %  'DrawNodes' - flag determining whether to draw nodes,
            %    default is true
            %
            %  'DrawBodies' - flag determining whether to draw bodies,
            %    default is true.
            %
            %  'OnlyNodes' - vector of indices of nodes to plot, only these
            %    nodes will be plotted. By default all nodes are plotted if
            %    'DrawNodes' is true. This setting is ignored if
            %    'DrawNodes' is false.
            %
            %  'Light' - flag determining whether to light the scene to
            %    give a more realistic look. Best used with 'solid' option.
            %
            %  'NodeSize' - scalar value indicating how large the nodes
            %    should be drawn (in the same units as the plot). If not
            %    supplied, the nodes are drawn with a size based on the
            %    axes limits.
            %
            %  'View' - viewpoint of plot in the same format as any of the
            %    single input styles of the 'view' function. See output of
            %    'help view' for more information. The contents of the
            %    'View' option is passed directly to the view function to
            %    set the axes view.
            %
            %  'FigPositionAndSize' - figure position and size as a four
            %    element vector in the same format as a figure's 'Position'
            %    property, i.e. [ x, y, w, h ] where x and y are the
            %    coordinates of the lower left corner of the figure and w
            %    and h are the height and width respectively.
            %
            % Output
            %
            %  hfig - handle to figure created
            %
            %  hax - handle to plot axes created
            %
            %
            
            if ~self.simComplete
                error ('Simulation data is not available for plotting')
            end
            
            options.PlotAxes = [];
            options.DrawLabels = false;
            options.AxLims = [];
            options.NodeSize = [];
            options.DrawMode = 'wireghost';
            options.DrawNodes = true;
            options.DrawBodies = true;
            options.Light = false;
            options.Title = true;
            options.OnlyNodes = [];
            options.ForceRedraw = false;
            options.DrawWaves = true;
            options.View = [];
            options.FigPositionAndSize = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkAxes (options.PlotAxes);
            
            if isempty (self.mBDynPostProc)
                self.mBDynPostProc = mbdyn.postproc (self.mBDynInputFile, self.mBDynSystem);
            end
            
            if isempty (options.OnlyNodes)
                options.OnlyNodes = self.mBDynPostProc.nNodes;
            end
            
            if options.DrawWaves
                
                wavedrawfcn = @(hax, tind) wavePlot3D ( self, ...
                                                        [], ...
                                                        'PlotAxes', hax, ...
                                                        'DomainSize', 'fromaxes', ...
                                                        'TimeIndex', tind );
                
                plotdata = self.mBDynPostProc.drawStep ( tind, ...
                                              'PlotAxes', self.drawAxesH, ...
                                              'DrawLabels', options.DrawLabels, ...
                                              'AxLims', options.AxLims, ...
                                              'NodeSize', options.NodeSize, ...
                                              'DrawMode', options.DrawMode, ...
                                              'DrawNodes', options.DrawNodes, ...
                                              'DrawBodies', options.DrawBodies, ...
                                              'Light', options.Light, ...
                                              'Title', options.Title, ...
                                              'OnlyNodes', options.OnlyNodes, ...
                                              'ForceRedraw', options.ForceRedraw, ...
                                              'ExternalDrawFcn', wavedrawfcn, ...
                                              'View', options.View, ...
                                              'FigPositionAndSize', options.FigPositionAndSize );
                                      
            else
                
                plotdata = self.mBDynPostProc.drawStep ( tind, ...
                                              'PlotAxes', self.drawAxesH, ...
                                              'DrawLabels', options.DrawLabels, ...
                                              'AxLims', options.AxLims, ...
                                              'NodeSize', options.NodeSize, ...
                                              'DrawMode', options.DrawMode, ...
                                              'DrawNodes', options.DrawNodes, ...
                                              'DrawBodies', options.DrawBodies, ...
                                              'Light', options.Light, ...
                                              'Title', options.Title, ...
                                              'OnlyNodes', options.OnlyNodes, ...
                                              'ForceRedraw', options.ForceRedraw, ...
                                              'View', options.View, ...
                                              'FigPositionAndSize', options.FigPositionAndSize  );
                                          
            end
            
        end
        
        
        function animate (self, varargin)
            % animate the wave energy system
            %
            % Syntax
            %
            % animate (wsobj, 'Parameter', value) 
            %
            % Input
            %
            %  wsobj - wsim.wecSim object
            %
            % Addtional optional arguments are supplied as parameter-value
            % pairs. The available options are:
            %
            %  'PlotAxes' - handle to plot axes in which to draw the
            %    trajectories. If not supplied a new figure and axes are
            %    created for the plot.
            %
            %  'AxLims' - (3 x 2) matrix containing limits the plot axes
            %    should be set to, each row is the x, y and z axes limits
            %    respectively. Uselful to focus on a particular region of
            %    the system. If not supplied, suitable axes limits will be
            %    determined based on the motion data.
            %
            %  'Title' - flag determining whether to add a title to the
            %    system plot. Default is true.
            %
            %  'DrawMode' - string indicating the drawing mode. Can be one
            %    of: 'solid', 'ghost', 'wire', 'wireghost'.
            %
            %  'DrawLabels' - flag indicating whether to draw node labels,
            %    default is false.
            %
            %  'DrawNodes' - flag determining whether to draw nodes,
            %    default is true
            %
            %  'DrawBodies' - flag determining whether to draw bodies,
            %    default is true.
            %
            %  'OnlyNodes' - vector of indices of nodes to plot, only these
            %    nodes will be plotted. By default all nodes are plotted if
            %    'DrawNodes' is true. This setting is ignored if
            %    'DrawNodes' is false.
            %
            %  'Light' - flag determining whether to light the scene to
            %    give a more realistic look. Best used with 'solid' option.
            %
            %  'NodeSize' - scalar value indicating how large the nodes
            %    should be drawn (in the same units as the plot). If not
            %    supplied, the nodes are drawn with a size based on the
            %    axes limits.
            %
            %  'ExternalDrawFcn' - function handle or string to function
            %    which  will be called after drawing the scene is complete.
            %    Intended to be used by external programs to add their own
            %    features to the scene. The funciton must take two
            %    arguments, the first is the handle to the axes contining
            %    the plot, the second is the time step index of the
            %    simulation.
            %
            %  'View' - viewpoint of plot in the same format as any of the
            %    single input styles of the 'view' function. See output of
            %    'help view' for more information. The contents of the
            %    'View' option is passed directly to the view function to
            %    set the axes view.
            %
            %  'FigPositionAndSize' - figure position and size as a four
            %    element vector in the same format as a figure's 'Position'
            %    property, i.e. [ x, y, w, h ] where x and y are the
            %    coordinates of the lower left corner of the figure and w
            %    and h are the height and width respectively.
            %
            
            if ~self.simComplete
                error ('No simulation results are available for plotting')
            end
            
            options.PlotAxes = [];
            options.DrawLabels = false;
            options.AxLims = [];
            options.PlotTrajectories = false;
            options.NodeSize = [];
            options.DrawMode = 'wireghost';
            options.DrawNodes = true;
            options.DrawBodies = true;
            options.Skip = 1;
            options.Light = false;
            options.VideoFile = [];
            options.VideoSpeed = 1;
            options.OnlyNodes = [];
            options.ExternalDrawFcn = [];
            options.View = [];
            options.DrawWaves = true;
            options.FigPositionAndSize  = [];
            
            options = parse_pv_pairs (options, varargin);
            
            self.checkAxes (options.PlotAxes);

            if isempty (self.mBDynPostProc)
                self.mBDynPostProc = mbdyn.postproc (self.mBDynInputFile, self.mBDynSystem);
            end

            if isempty (options.OnlyNodes)
                options.OnlyNodes = 1:self.mBDynPostProc.nNodes;
            end

            if options.DrawWaves

                wavedrawfcn = @(hax, tind) wavePlot3D ( self, ...
                                                        [], ...
                                                        'PlotAxes', hax, ...
                                                        'DomainSize', 'fromaxes', ...
                                                        'TimeIndex', tind );

                self.mBDynPostProc.animate ( ...
                                              'PlotAxes', self.drawAxesH, ...
                                              'DrawLabels', options.DrawLabels, ...
                                              'AxLims', options.AxLims, ...
                                              'NodeSize', options.NodeSize, ...
                                              'DrawMode', options.DrawMode, ...
                                              'DrawNodes', options.DrawNodes, ...
                                              'DrawBodies', options.DrawBodies, ...
                                              'Light', options.Light, ...
                                              'VideoFile', options.VideoFile, ...
                                              'VideoSpeed', options.VideoSpeed, ...
                                              ...'Title', options.Title, ...
                                              'OnlyNodes', options.OnlyNodes, ...
                                              ...'ForceRedraw', options.ForceRedraw, ...
                                              'ExternalDrawFcn', wavedrawfcn, ...
                                              'View', options.View, ...
                                              'FigPositionAndSize', options.FigPositionAndSize );

            else

                self.mBDynPostProc.animate (  ...
                                              'PlotAxes', self.drawAxesH, ...
                                              'DrawLabels', options.DrawLabels, ...
                                              'AxLims', options.AxLims, ...
                                              'NodeSize', options.NodeSize, ...
                                              'DrawMode', options.DrawMode, ...
                                              'DrawNodes', options.DrawNodes, ...
                                              'DrawBodies', options.DrawBodies, ...
                                              'Light', options.Light, ...
                                              'VideoFile', options.VideoFile, ...
                                              'VideoSpeed', options.VideoSpeed, ...
                                              ...'Title', options.Title, ...
                                              'OnlyNodes', options.OnlyNodes, ...
                                              ...'ForceRedraw', options.ForceRedraw, ...
                                              'View', options.View, ...
                                              'FigPositionAndSize', options.FigPositionAndSize  );

            end

        end
        
        
        function [hsurf, hax] = wavePlot3D (self, t, varargin)
            % plot the wave elevation as a 3D surface
            
            options.PlotAxes = [];
            options.NumPointsX = self.wavePlotNumPointsX;
            options.NumPointsY = self.wavePlotNumPointsY;
            options.DomainSize = self.wavePlotDomainSize;
            options.DomainCorner = [];
            options.TimeIndex = [];
            options.ForceRedraw = false;
            
            options = parse_pv_pairs (options, varargin);
            
            check.isLogicalScalar (options.ForceRedraw, true, 'ForceRedraw');
            
            self.checkAxes (options.PlotAxes);
            
            hax = self.drawAxesH;
            
            if ~isempty (options.TimeIndex)
                % override the t value with logged time from index
                check.isPositiveScalarInteger (options.TimeIndex, true, 'TimeIndex', false);
                t = self.logger.data.Time (options.TimeIndex);
            end
            
%             xlim = get (hax, 'Xlim');
%             ylim = get (hax, 'Ylim');
            
            if ischar (options.DomainSize)
                switch options.DomainSize
                    
                    case 'fromaxes'
                        
                        xlim = get (hax, 'Xlim');
                        ylim = get (hax, 'Ylim');
                        options.DomainSize = [ xlim(2) - xlim(1), ylim(2) - ylim(1) ];
                        options.DomainCorner = [ xlim(1), ylim(1) ];
                        
                end
            end
            
            % the following test is supposed to capture when we are doing a
            % brand new plot, or changing the domain of the existing plot
            if isempty (self.wavePlotDomainSize) && isempty (options.DomainSize)
                
                options.DomainSize = [ 10, 10 ];
                
            elseif isscalar (options.DomainSize)
                
                options.DomainSize = [options.DomainSize, options.DomainSize];
                
            end
            
            if isempty (options.DomainCorner)
                options.DomainCorner = [ -options.DomainSize(1)/2, -options.DomainSize(2)/2 ];
            end
            
            % check if the desired plot parameters have changed at all,
            % requiring the surface to be re-plotted
            if isempty (self.waveSurfH) ...
                || isempty (self.wavePlotDomainSize) ...
                || options.NumPointsX ~= self.wavePlotNumPointsX ...
                || options.NumPointsY ~= self.wavePlotNumPointsY ...
                || ~all (options.DomainSize == self.wavePlotDomainSize) ...
                || options.ForceRedraw ...
                || ~ishghandle (self.waveSurfH)
                
                self.wavePlotNumPointsX = options.NumPointsX;
                self.wavePlotNumPointsY = options.NumPointsY;
                self.wavePlotDomainSize = options.DomainSize;
                
                % make the grid
                x = linspace ( options.DomainCorner(1), ...
                               options.DomainCorner(1) + options.DomainSize(1), ...
                               self.wavePlotNumPointsX );
                           
                y = linspace ( options.DomainCorner(2), ...
                               options.DomainCorner(2) + options.DomainSize(2), ...
                               self.wavePlotNumPointsY );
                
                [self.wavePlotX, self.wavePlotY] = meshgrid (x, y);
                
                wavecolour = [0.3010, 0.7450, 0.9330];
                
                if ~isempty (self.waveSurfH) && ishghandle (self.waveSurfH)
                    delete (self.waveSurfH);
                end

                hold on
                self.waveSurfH = surf ( hax, ...
                                        self.wavePlotX, ...
                                        self.wavePlotY, ...
                                        self.hydroSystem.waves.waveElevationGrid (t, self.wavePlotX, self.wavePlotY), ...
                                        'EdgeColor', wavecolour, ...
                                        'FaceColor', wavecolour, ...
                                        'EdgeAlpha', 0.25, ...
                                        'FaceAlpha', 0.25 );
                hold off
                                    
            else
                % there is an existing plot
                
                % recalculate Z
                Z = self.hydroSystem.waves.waveElevationGrid (t, self.wavePlotX, self.wavePlotY);
            
                % update the z data for the surface plot
                set ( self.waveSurfH, 'Zdata', Z );
                
                % drawnow () % don't know if we need a drawnow call here 
            end
            
            if nargout > 0
                hsurf = self.waveSurfH;
            end
            
        end
        
    end
    
    methods
        % getter/setter methods go here
        
        function set.loggingSettings (self, newsettings)
            
            assert (isa (newsettings, 'wsim.loggingSettings'), ...
                'loggingSettings property must be an wsim.loggingSettings object');
            
            self.loggingSettings = newsettings;
            
        end
        
        function set.wavePlotNumPointsX (self, n)
            
            check.isPositiveScalarInteger (n, true, 'wavePlotNumPointsX', true);
            
            self.wavePlotNumPointsX = n;
            
        end
        
        function set.wavePlotNumPointsY (self, n)
            
            check.isPositiveScalarInteger (n, true, 'wavePlotNumPointsY', true);
            
            self.wavePlotNumPointsY = n;
            
        end
        
    end
    
    methods (Access = private)
        
        function [pos, vel, accel] = getMotion (self, mb)
            
            eul = zeros (3, self.nMBDynNodes);
            
            R = mb.GetRot();
            for Rind = 1:size (R,3)
%                 om = mbdyn.pre.orientmat ('orientation', R(:,:,Rind));
%                 eul(1:3,Rind) = om.euler123 ();

                eul(1:3,Rind) = self.euler123(R(:,:,Rind));

            end
            
            pos = [ mb.NodePositions();
                    eul ];
            
%             pos = [ mb.NodePositions();
%                     mb.GetRot() ];

            vel = [ mb.NodeVelocities();
                    mb.NodeOmegas() ];
            
            accel = [ mb.NodeAccelerations();
                      mb.NodeAngularAccels() ];
            
        end
        
        function eul = euler123 (self, R)
            % returns the extrinsic euler123 angles corresponding to the
            % orientation matrix

            alpha = -atan2 ( R(2,3) , R(3,3) );

            beta = atan2 ( R(1,3) ...
                            , ( cos (alpha)*R(3,3) - sin (alpha)*R(2,3) ) );

            gamma = atan2 ( ( cos (alpha)*R(2,1) - sin (alpha)*R(3,1) ) ...
                           , ( cos (alpha)*R(2,2) - sin (alpha)*R(3,2) ) );

            eul = [ alpha;
                    beta;
                    gamma ];
                
        end
        
        function forces_and_moments = applyPTOForces (self, forces_and_moments)
            
            % get the PTO forces and moments and put them in the corect
            % places in the forces_and_moments matrix, based on an index
            % mapping done by the method mapPTOForceInds(), which is called
            % by the prepare method before simulation can be started.
            for ptoind = 1:size (self.ptoIndexMap, 1)
                
                % Each PTO has a reference node and another node, the PTO
                % forces/moments are applied between the two nodes. The
                % returned PTO force is applied to the reference node, and
                % the opposite force is applied to the other node
                
                ptoForceAndTorque = self.powerTakeOffs{ptoind}.forceAndMoment (self.lastTime);

                forces_and_moments (:,self.ptoIndexMap(ptoind,1)) = ...
                    forces_and_moments (:,self.ptoIndexMap(ptoind,1)) + ptoForceAndTorque(:,1);
                
                forces_and_moments (:,self.ptoIndexMap(ptoind,2)) = ...
                    forces_and_moments (:,self.ptoIndexMap(ptoind,2)) + ptoForceAndTorque(:,2);
                
            end
            
            
        end
        
        function startPTOs (self, siminfo)
            
            for ptoind = 1:size (self.ptoIndexMap, 1)
                
                self.powerTakeOffs{ptoind}.start (siminfo);
                
            end 
            
        end
        

        function forces_and_moments = applyHydroForces (self, forces_and_moments, time, pos, vel, accel)
                    
            if self.hydroMotionSyncStepCount == self.hydroMotionSyncStep
                
                % calculate hydrodynamic forces based on the motion
                [hydroforces, out] = self.hydroSystem.hydroForces (time, ...
                                                pos(:,self.hydroNodeIndexMap), ...
                                                vel(:,self.hydroNodeIndexMap), ...
                                                accel(:,self.hydroNodeIndexMap) );

                

                % store the last generated forces fo later use (e.g. logging)
                self.lastForceHydro = hydroforces;
                self.lastForceExcitation = out.F_Excit;
                self.lastForceExcitationRamp = out.F_ExcitRamp;
                self.lastForceExcitationLin = out.F_ExcitLin;
                self.lastForceExcitationNonLin = out.F_ExcitNonLin;
                self.lastForceRadiationDamping = out.F_RadiationDamping;
                self.lastForceRestoring = out.F_Restoring;
                self.lastForceMorrison = out.F_MorrisonElement;
                self.lastForceViscousDamping = out.F_ViscousDamping;
                self.lastForceAddedMassUncorrected = out.F_AddedMass;
                
                % reset step count
                self.hydroMotionSyncStepCount = 1;
            
            else
                
                % just return the last set of hydro forces calculated
                hydroforces = self.lastForceHydro;
                
                self.hydroMotionSyncStepCount = self.hydroMotionSyncStepCount + 1;
                
            end
            
            % add hydrodynamic forces to the correct nodes (which are the
            % nodes attached to bodies with hydrodynamic interaction). This
            % is done using a mapping of the indexes created by the
            % mapHydroForceInds method, which is called by the 'prepare'
            % method
            for ind = 1:size (self.hydroNodeIndexMap, 1)
                forces_and_moments (:,self.hydroNodeIndexMap(ind,1)) = ...
                    forces_and_moments (:,self.hydroNodeIndexMap(ind,1)) + hydroforces(:,ind);
            end
            
        end
        
        function mapPTOForceInds (self)
            % gets the node indices for application of power take-off forces
            %
            % Syntax
            %
            % wsim.mapPTOForceInds (ws)
            %
            % Description
            %
            % mapPTOForceInds creates a matrix whose values give the
            % appropriate indices of the force and moment matrix for the
            % whole system to which to add the forces from each PTO. Not
            % every node in the system is necessarily accessible via the
            % interface, only those associated with the MBDyn external
            % structural forces element, so this is actually a mapping of
            % the index of the structural external force nodes to the PTO
            % nodes. The mapping is done by comparing the mbdyn structural
            % nodes associated with the PTO objects with all the nodes in
            % the external sturctural force element.
            %
            % mapPTOForceInds updates the ptoIndexMap of the wsim.wecSim
            % object to contain a two column matrix. Each row corresponds
            % to each of the PTO objects. The first column is the index of
            % the column of the force_and_moment matrix corresponding to
            % the reference node of each PTO, the second column is the
            % index of the column of the force_and_moment matrix
            % corresponding to the other node of each PTO.
            %
            % Input
            %
            %  ws - wsim.wecSim object
            %
            %
            
            strinfo = self.mBDynSystem.externalStructuralInfo ();
            
            % preallocate the index map for all the nodes
            self.ptoIndexMap = zeros (numel (self.powerTakeOffs), 2);
            
            for ptoind = 1:numel(self.powerTakeOffs)
                % loop through the power take-off objects, for each one
                % match up the PTO nodes with the nodes in the whole
                % system. Each PTO has a reference node and another node,
                % the PTO force/moment are applied between the two nodes.
                % The returned PTO force is applied to the reference node,
                % and the opposite force is applied to the other node
                %
                
                foundptonode = false;
                for nodeind = 1:numel (strinfo.Nodes)
                    
                    if strinfo.Nodes{nodeind} == self.powerTakeOffs{ptoind}.referenceNode
                        
                        foundptonode = true;
                        
                        self.ptoIndexMap(ptoind, 1) = nodeind;
                        
                        break;
                    
                    end
                
                end
                
                if foundptonode == false
                    error ('Could not find reference node for PTO %d in MBDyn system structural external force element', ptoind);
                end
                
                foundptonode = false;
                for nodeind = 1:numel (strinfo.Nodes)
                    
                    if strinfo.Nodes{nodeind} == self.powerTakeOffs{ptoind}.otherNode
                        
                        foundptonode = true;
                        
                        self.ptoIndexMap(ptoind, 2) = nodeind;
                    
                        break;
                        
                    end
                
                end
                
                if foundptonode == false
                    error ('Could not find non-reference node for PTO %d in MBDyn system structural external force element', ptoind);
                end
                
            end
            
        end
        
        function mapHydroForceInds (self)
            % gets the node indices for application of hydrodynamic forces
            %
            % Syntax
            %
            % wsim.mapHydroForceInds (ws)
            %
            % Description
            %
            % mapPTOForceInds creates a matrix whose values give the
            % appropriate indices of the force and moment matrix for the
            % whole system to which to add the forces due to hydrodynamic
            % interaction. Not every node in the system is necessarily
            % accessible via the interface, only those associated with the
            % MBDyn external structural forces element, so this is actually
            % a mapping of the index of the structural external force nodes
            % to the hydro body nodes. The mapping is done by comparing the
            % mbdyn structural nodes associated with the hydrobody objects
            % with all the nodes in the external structural force element.
            %
            % mapHydroForceInds updates the hydroNodeIndexMap of the
            % wsim.wecSim object to contain a column vector containing the
            % index of the column of the force_and_moment matrix
            % corresponding to each hydrobody in the system.
            %
            % Input
            %
            %  ws - wsim.wecSim object
            %
            
            
            strinfo = self.mBDynSystem.externalStructuralInfo ();
            
            self.hydroNodeIndexMap = zeros (numel (self.hydroSystem.bodyMBDynNodes), 1);
            
            for hydronodeind = 1:numel(self.hydroSystem.bodyMBDynNodes)
                
                foundnode = false;
                for nodeind = 1:numel (strinfo.Nodes)
                    
                    if strinfo.Nodes{nodeind} == self.hydroSystem.bodyMBDynNodes{hydronodeind}
                        
                        foundnode = true;
                        
                        self.hydroNodeIndexMap(hydronodeind, 1) = nodeind;
                        
                        break;
                    
                    end
                
                end
                
                if foundnode == false
                    error ('Could not find corresponding external structural element node for hydro node %d', hydronodeind);
                end
                
            end
            
        end
        

        function initDataStructures (self)
            % intialise logging data structures

            extforceinfo = self.mBDynSystem.externalStructuralInfo ();

            % make a clean logger object to store simulation data
            self.logger = wsim.logger ();
            
            % preallocate data vectors for the predicted number of time
            % steps
            nsteps = (self.hydroSystem.simu.endTime - self.hydroSystem.simu.startTime) ./ self.hydroSystem.simu.dt + 1;
            
            % always store time, as is is independent variabl for other
            % vars in logger object
            self.logger.addVariable ( 'Time', [1, 1], ...
                                      'Desc', 'main hydrodynamic/multibody time step', ...
                                      'AxisLabel', 'Time [s]', ...
                                      'Pre', nsteps );
            
            if self.loggingSettings.positions
                self.logger.addVariable ( 'Positions', [3, extforceinfo.NNodes], ...
                                          'Desc', 'cartesian positions of all structural external nodes', ...
                                          'AxisLabel', 'External Struct Nodes Positions [m]', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.angularPositions
                self.logger.addVariable ( 'AngularPositions', [3, extforceinfo.NNodes], ...
                                          'Desc', 'angular positions of all structural external nodes', ...
                                          'AxisLabel', 'External Struct Nodes Euler Angles [rad]', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.velocities
                self.logger.addVariable ( 'Velocities', [3, extforceinfo.NNodes], ...
                                          'Desc', 'cartesian velocities of all structural external nodes', ...
                                          'AxisLabel', 'External Struct Nodes Velocities [ms^{-1}]', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.angularVelocities
                self.logger.addVariable ( 'AngularVelocities', [3, extforceinfo.NNodes], ...
                                          'Desc', 'angular velocities of all structural external nodes', ...
                                          'AxisLabel', 'External Struct Nodes Angular Velocities [rads^{-1}]', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.accelerations
                self.logger.addVariable ( 'Accelerations', [3, extforceinfo.NNodes], ...
                                          'Desc', 'cartesian accelerations of all structural external nodes', ...
                                          'AxisLabel', 'External Struct Nodes Accelerations [ms^{-2}]', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.angularAccelerations
                self.logger.addVariable ( 'AngularAccelerations', [3, extforceinfo.NNodes], ...
                                          'Desc', 'angular accelerations of all structural external nodes', ...
                                          'AxisLabel', 'External Struct Nodes Angular Accelerations [rad s^{-2}]', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.nodeForcesAndMoments || self.loggingSettings.forceAddedMass
                self.logger.addVariable ( 'NodeForcesAndMoments', [6, extforceinfo.NNodes], ...
                                          'Desc', 'sum of all forces for all nodes with external forces', ...
                                          'AxisLabel', 'Total Forces [N] and Moments [Nm] on External Struct Nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.nodeForcesAndMomentsUncorrected || self.loggingSettings.forceAddedMass
                self.logger.addVariable ( 'NodeForcesAndMomentsUncorrected', [6, extforceinfo.NNodes], ...
                                          'Desc', 'sum of all forces (with uncorrected added mass forces) for all external structural nodes with external forces', ...
                                          'AxisLabel', 'Uncorrected Added Mass Forces [N] and Moments [Nm] on Nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            
            if self.loggingSettings.forceHydro
                self.logger.addVariable ( 'ForceHydro', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'sum of all hydrodynamic forces for all hydro nodes', ...
                                          'AxisLabel', 'Total Hydro Forces [N] and Moments [Nm] on Nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceExcitation
                self.logger.addVariable ( 'ForceExcitation', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'sum of linear and nonlinear hydrodynamic excitation force for all hydro nodes', ...
                                          'AxisLabel', 'Total Hydro Excitation Forces [N] and Moments [Nm] on Nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceExcitationRamp
                self.logger.addVariable ( 'ForceExcitationRamp', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'sum of linear and nonlinear hydrodynamic excitation force for all hydro nodes, but with a ramp function applied', ...
                                          'AxisLabel', 'Ramped Hydro Excitation Forces [N] and Moments [Nm] on Nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceExcitationLin
                self.logger.addVariable ( 'ForceExcitationLin', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'linear hydrodynamic excitation forces for all hydro nodes', ...
                                          'AxisLabel', 'Linear Hydro Excitation Forces [N] and Moments [Nm] on Nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceExcitationNonLin
                self.logger.addVariable ( 'ForceExcitationNonLin', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'nonlinear hydrodynamic excitation forces for all hydro nodes', ...
                                          'AxisLabel', 'Nonlinear Hydro Excitation Forces [N] and Moments [Nm] on Nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceRadiationDamping
                self.logger.addVariable ( 'ForceRadiationDamping', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'hydrodynamic radiation and damping forces for all hydro nodes', ...
                                          'AxisLabel', 'Hydro Radiation Forces [N] and Moments [Nm] on Nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceRestoring
                self.logger.addVariable ( 'ForceRestoring', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'hydrodynamic restoring forces for all hydro nodes', ...
                                          'AxisLabel', 'Hydro Restoring Forces [N] and Moments [Nm] on Nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceMorrison
                self.logger.addVariable ( 'ForceMorrison', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'morrison forces for all hydro nodes', ...
                                          'AxisLabel', 'Morrison Forces [N] and Moments [Nm] on Nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceViscousDamping
                self.logger.addVariable ( 'ForceViscousDamping', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'hydrodynamic viscous damping forces for all hydro nodes', ...
                                          'AxisLabel', 'Viscous Damping Forces [N] and Moments [Nm] on Nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceAddedMass
                self.logger.addVariable ( 'ForceAddedMass', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'hydrodynamic added mass forces for all hydro nodes', ...
                                          'AxisLabel', 'Added Mass Forces [N] and Moments [Nm] on Nodes ', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceAddedMassUncorrected
                self.logger.addVariable ( 'ForceAddedMassUncorrected', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'uncorrected hydrodynamic added mass forces for all hydro nodes', ...
                                          'AxisLabel', 'Uncorrected Added Mass Forces [N] and Moments [Nm] on Nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.powerTakeOffInternal
                
                for ptoind = 1:numel (self.powerTakeOffs)
                    
                    self.powerTakeOffs{ptoind}.loggingSetup (self.logger);
                    
                end
                
            end

        end
        
        
        function logData (self)
            % log data after a time step has completed

            % always store time, as is is independent variable for other
            % vars in logger object
            self.logger.logVal ( 'Time', self.lastTime );
            
            if self.loggingSettings.positions
                self.logger.logVal ( 'Positions', self.lastPositions );
            end
            
            if self.loggingSettings.angularPositions
                self.logger.logVal ( 'AngularPositions', self.lastAngularPositions );
            end
            
            if self.loggingSettings.velocities
                self.logger.logVal ( 'Velocities', self.lastVelocities );
            end
            
            if self.loggingSettings.angularVelocities
                self.logger.logVal ( 'AngularVelocities', self.lastAngularVelocities );
            end
            
            if self.loggingSettings.angularAccelerations
                self.logger.logVal ( 'AngularAccelerations', self.lastAngularAccelerations );
            end
            
            if self.loggingSettings.accelerations
                self.logger.logVal ( 'Accelerations', self.lastAccelerations );
            end
            
            if self.loggingSettings.nodeForcesAndMomentsUncorrected ...
                    || self.loggingSettings.forceAddedMass
                self.logger.logVal ( 'NodeForcesAndMomentsUncorrected', self.lastNodeForcesAndMomentsUncorrected );
            end
            
            if self.loggingSettings.forceHydro
                self.logger.logVal ( 'ForceHydro', self.lastForceHydro );
            end
            
            if self.loggingSettings.forceExcitation
                self.logger.logVal ( 'ForceExcitation', self.lastForceExcitation );
            end
            
            if self.loggingSettings.forceExcitationRamp
                self.logger.logVal ( 'ForceExcitationRamp', self.lastForceExcitationRamp );
            end
            
            if self.loggingSettings.forceExcitationLin
                self.logger.logVal ( 'ForceExcitationLin', self.lastForceExcitationLin );
            end
            
            if self.loggingSettings.forceExcitationNonLin
                self.logger.logVal ( 'ForceExcitationNonLin', self.lastForceExcitationNonLin );
            end
            
            if self.loggingSettings.forceRadiationDamping
                self.logger.logVal ( 'ForceRadiationDamping', self.lastForceRadiationDamping );
            end
            
            if self.loggingSettings.forceRestoring
                self.logger.logVal ( 'ForceRestoring', self.lastForceRestoring );
            end
            
            if self.loggingSettings.forceMorrison
                self.logger.logVal ( 'ForceMorrison', self.lastForceMorrison );
            end
            
            if self.loggingSettings.forceViscousDamping
                self.logger.logVal ( 'ForceViscousDamping', self.lastForceViscousDamping );
            end
            
            if self.loggingSettings.forceAddedMass || self.loggingSettings.forceAddedMassUncorrected
                self.logger.logVal ( 'ForceAddedMassUncorrected', self.lastForceAddedMassUncorrected );
            end
            
            
        end
        
        function advanceStep (self)
            
            self.logData ();
            
            % accept the last data into the time history of solutions
            % for the hydrodynamic system and advance
            self.hydroSystem.advanceStep ( self.lastTime, ...
                                           [ self.lastVelocities(:,self.hydroNodeIndexMap(:,1)); ...
                                             self.lastAngularVelocities(:,self.hydroNodeIndexMap(:,1)) ], ...
                                           [ self.lastAccelerations(:,self.hydroNodeIndexMap(:,1)); ...
                                             self.lastAngularAccelerations(:,self.hydroNodeIndexMap(:,1)) ] );
            
            for ptoind = 1:numel (self.powerTakeOffs)
                % PTOs handle their own data logging which is triggered by
                % calling advanceStep on each PTO object
                self.powerTakeOffs{ptoind}.advanceStep (self.lastTime);
                    
            end
            
        end
        
        function finishPTOs (self)
            % tell all the PTOs that the simulation is complete
            
            for ptoind = 1:numel (self.powerTakeOffs)

                self.powerTakeOffs{ptoind}.finish (self.lastTime);
                    
            end
            
        end
        
        function displayLastNLinesOfFile (self, filename, nlines)
            % display the last n lines of a text file on the command line
            
            lines = cell (nlines,1);
            
            fid = fopen (filename);
            
            if fid ~= -1
                
                CC = onCleanup (@() fclose (fid));

                line1ind = nlines + 1;
                
                while ~feof(fid)
                    
                    line1ind = max (1, line1ind - 1);
                    
                    lines = circshift (lines, -1, 1);
                    
                    lines{end,1} = fgetl(fid);
                    
                end
                
                % only print the lines we actually read in if there's less
                % than the max allowed
                for ind = line1ind:nlines
                    fprintf (1, '%s\n', lines{ind});
                end
            
            end
            
        end
        
        function makeAxes (self)
            % create axes and transform object
            
            figure;
            self.drawAxesH = axes;
            
        end
        
        function checkAxes (self, hax)
            % checks if there is a valid set of axes to plot to, and if
            % not, creates them
            
            % try to figure out if there is a valid set of axes to plot to
            if isempty (hax)
                
                % plot in the existing axes if possible as no new axes
                % supplied
                if mbdyn.pre.base.isAxesHandle (self.drawAxesH)
                    % use existing axes
                    
                    if ~ishghandle (self.drawAxesH)
                        % axes no longer exists, or figure has been closed
                        self.drawAxesH = [];
%                         self.deleteAllDrawnObjects ();
                        % make a new set of axes to draw to
                        self.makeAxes ();
                    end
                    
                elseif isempty (self.drawAxesH)
                    % make figure and axes
                    self.makeAxes ();
%                     self.needsRedraw = true;
                    
                else
                    error ('drawAxesH property is not empty or an axes handle');
                end
            
            elseif mbdyn.pre.base.isAxesHandle (hax)
                % an axes has been supplied, so we plot to this new axes
                
                if isoctave
                    drawnow ();
                end
                
                if ~ishghandle (hax)
                    error ('provided axes object is not valid');
                end
                self.drawAxesH = hax;
               
%                 self.needsRedraw = true;
                
            else
                
                error ('hax must be an axes handle or empty');
                
            end
            
        end
        
    end
    
end