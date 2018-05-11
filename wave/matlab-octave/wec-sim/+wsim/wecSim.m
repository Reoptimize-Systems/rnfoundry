classdef wecSim < handle
    
    properties (GetAccess = public, SetAccess = private)
        
        powerTakeOffs;
        forceRelTol;
        forceAbsTol;
        mBDynOutputFile;
        mBDynInputFile;

        loggingSettings;
        
    end
    
    properties (GetAccess = private, SetAccess = private)
        
        ptoIndexMap;
        hydroNodeIndexMap;
        mBDynSystem;
        hydroSystem;
        readyToRun;
        hydroMotionSyncStep = 1;
        hydroMotionSyncStepCount = 1;
        
        
        logger;
        nMBDynNodes;
        
        % logging variables
        lastTime;                   % last value of simulation time
        lastPositions;              % last calculated positions and euler angles
        lastVelocities;             % last calculated velocites and angular velocities
        lastAccelerations;          % last calculated accelerations and angular accelerations
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
        
    end
    
    methods
        
        function self = wecSim (hsys, mbsys, varargin)
            % wsim.wecSim constructor
            %
            
            options.PTOs = {};
%             options.StartTime = [];
%             options.EndTime = [];
            options.NEMOHSim = [];
            options.MBDynInputFile = '';
            options.LoggingSettings = wsim.loggingSettings ();
            options.MBDynVerbosity = 0;
            
            options = parse_pv_pairs (options, varargin);
            
            assert (isa (hsys, 'wsim.hydroSystem'), ...
                'hsys must be an wsim.hydroSystem object');
            
            assert (isa (mbsys, 'mbdyn.pre.system'), ...
                'mbsys must be an mbdyn.pre.system object');
            
            if ~isempty (options.NEMOHSim)
                assert (isa (options.NEMOHSim, 'nemoh.simulation'), ...
                    'NEMOHSim must be a nemoh.simulation object' ); 
            end
            
            if isempty (options.MBDynInputFile)
                options.MBDynInputFile = fullfile (hsys.simu.caseDir, 'mbdyn_input_file.mbd');
            else
                assert (ischar (options.MBDynInputFile), ...
                    'MBDynInputFile must be string containing the file path where the MBDyn input file should be generated');
            end
            
            self.loggingSettings = options.LoggingSettings;
            self.mBDynSystem = mbsys;
            self.hydroSystem = hsys;
            self.readyToRun = false;
            self.mBDynInputFile = options.MBDynInputFile;
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
            
            self.readyToRun = true;
            
        end
        
        function datalog = run (self, varargin)
            % Run the WEC simulation
            
            options.OutputFilePrefix = fullfile (self.hydroSystem.simu.caseDir, 'output', 'wsim');
            options.Verbosity = 0;
            options.AbsForceTolerance = 100;
            options.RelForceTolerance = 1e-5;
            options.MinIterations = 0;
            options.MaxIterations = self.mBDynSystem.problems{1}.maxIterations;
            options.TimeExecution = false;
            options.HydroMotionSyncSteps = 1;
            
            options = parse_pv_pairs (options, varargin);
            
            check.isPositiveScalarInteger (options.Verbosity, true, 'Verbosity', true);
            
            assert (ischar (options.OutputFilePrefix), ...
                'OutputFilePrefix must be a string');
            
            % check for positive numeric scalar
            check.isNumericScalar (options.AbsForceTolerance, true, 'AbsForceTolerance', 1);
            check.isNumericScalar (options.RelForceTolerance, true, 'RelForceTolerance', 1);
            
            check.isPositiveScalarInteger (options.MaxIterations, true, 'Verbosity');
            
            assert (options.MinIterations <= options.MaxIterations, ...
                'MinIterations (%d) is not less than or equal to MaxIterations (%d)', ...
                options.MinIterations,  options.MaxIterations)
            
            check.isLogicalScalar (options.TimeExecution, true, 'TimeExecution');
            
            if self.readyToRun == false
                error ('Simulation is not ready to run, have you run ''prepare'' yet?');
            end
            
            check.isPositiveScalarInteger (options.HydroMotionSyncSteps, true, 'HydroMotionSyncSteps')
            
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
            
            % start the PTO components
            self.startPTOs (siminfo);
            
            % generate input file and start mbdyn
            
            % clear out previous files
            delete ([options.OutputFilePrefix, '.*']);
            
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
            
            % ensure MBCNodal is destroyed in the event of a problem
            % (closes communication to MBDyn and tells it to quit so
            % sockets and so on are also cleaned up)
            CC = onCleanup (@() delete (mb));
            
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
            self.lastPositions = pos;
            self.lastVelocities = vel;
            self.lastAccelerations = accel;
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
                self.lastPositions = pos;
                self.lastVelocities = vel;
                self.lastAccelerations = accel;
                self.lastNodeForcesAndMomentsUncorrected = forces_and_moments;
                
                self.advanceStep ();
                
                ind = ind + 1;
                
            end
            
            fprintf ( 1, 'Reached time %f, in %d steps (of an intended %d), postprocessing ...\n', ...
                      self.lastTime, ...
                      ind-1, ...
                      self.logger.info.Time.PreallocatedLogLength );
            
            if self.loggingSettings.nodeForcesAndMoments || self.loggingSettings.forceAddedMass
                
                self.logger.setSeries ('NodeForcesAndMoments', self.logger.data.NodeForcesAndMomentsUncorrected);

                % the total forces on hydrodynamic bodies must be corrected
                [ self.logger.data.NodeForcesAndMoments(:,self.hydroNodeIndexMap(:,1),:), ...
                  F_added_mass ] = ...
                    correctAddedMassForce ( self.hydroSystem, ...
                                            self.logger.data.NodeForcesAndMomentsUncorrected(:,self.hydroNodeIndexMap(:,1),:), ...
                                            self.logger.data.ForceAddedMassUncorrected, ...
                                            self.logger.data.Accelerations );

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
            
        end
        
    end
    
    methods
        % getter/setter methods go here
        
        function set.loggingSettings (self, newsettings)
            
            assert (isa (newsettings, 'wsim.loggingSettings'), ...
                'loggingSettings property must be an wsim.loggingSettings object');
            
            self.loggingSettings = newsettings;
            
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
                                      'Pre', nsteps );
            
            if self.loggingSettings.positions
                self.logger.addVariable ( 'Positions', [6, extforceinfo.NNodes], ...
                                          'Desc', 'cartesian and angular positions of all structural external nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.velocities
                self.logger.addVariable ( 'Velocities', [6, extforceinfo.NNodes], ...
                                          'Desc', 'cartesian and angular velocities of all structural external nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.accelerations
                self.logger.addVariable ( 'Accelerations', [6, extforceinfo.NNodes], ...
                                          'Desc', 'cartesian and angular accelerations of all structural external nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.nodeForcesAndMoments || self.loggingSettings.forceAddedMass
                self.logger.addVariable ( 'NodeForcesAndMoments', [6, extforceinfo.NNodes], ...
                                          'Desc', 'sum of all forces for all nodes with external forces', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.nodeForcesAndMomentsUncorrected || self.loggingSettings.forceAddedMass
                self.logger.addVariable ( 'NodeForcesAndMomentsUncorrected', [6, extforceinfo.NNodes], ...
                                          'Desc', 'sum of all forces (with uncorrected added mass forces) for all nodes with external forces', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            
            if self.loggingSettings.forceHydro
                self.logger.addVariable ( 'ForceHydro', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'sum of all hydrodynamic forces for all hydro nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceExcitation
                self.logger.addVariable ( 'ForceExcitation', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'sum of linear and nonlinear hydrodynamic excitation force for all hydro nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceExcitationRamp
                self.logger.addVariable ( 'ForceExcitationRamp', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'sum of linear and nonlinear hydrodynamic excitation force for all hydro nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceExcitationLin
                self.logger.addVariable ( 'ForceExcitationLin', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'linear hydrodynamic excitation forces for all hydro nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceExcitationNonLin
                self.logger.addVariable ( 'ForceExcitationNonLin', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'nonlinear hydrodynamic excitation forces for all hydro nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceRadiationDamping
                self.logger.addVariable ( 'ForceRadiationDamping', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'hydrodynamic radiation and damping forces for all hydro nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceRestoring
                self.logger.addVariable ( 'ForceRestoring', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'hydrodynamic restoring forces for all hydro nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceMorrison
                self.logger.addVariable ( 'ForceMorrison', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'morrison forces for all hydro nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceViscousDamping
                self.logger.addVariable ( 'ForceViscousDamping', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'hydrodynamic viscous damping forces for all hydro nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceAddedMass
                self.logger.addVariable ( 'ForceAddedMass', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'hydrodynamic added mass forces for all hydro nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if self.loggingSettings.forceAddedMassUncorrected
                self.logger.addVariable ( 'ForceAddedMassUncorrected', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'uncorrected hydrodynamic added mass forces for all hydro nodes', ...
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
            
            if self.loggingSettings.velocities
                self.logger.logVal ( 'Velocities', self.lastVelocities );
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
                                           self.lastVelocities(:,self.hydroNodeIndexMap(:,1)), ...
                                           self.lastAccelerations(:,self.hydroNodeIndexMap(:,1)) );
            
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
        
    end
    
end