classdef wecsim < handle
    
    properties
        
%         hydroForceSystem;
        powerTakeOffs;
%         multibodySystem;
        
        forceRelTol;
        forceAbsTol;
        
    end
    
    properties (GetAccess = private, SetAccess = private)
        ptoIndexMap;
        hydroNodeIndexMap;
        mBDynSystem;
        hydroSystem;
        readyToRun;
        mBDynInputFile;
        logger;
    end
    
    methods
        
        function self = wecsim (hsys, mbsys, varargin)
            
            options.PTO = {};
            options.StartTime = [];
            options.EndTime = [];
            options.NEMOHSim = [];
            options.MBDynInputFile = '';
            
            options = parse_pv_pairs (options, varargin);
            
            assert (isa (hsys, 'wsim.hydrosys'), ...
                'hsys must be an wsim.hydrosys object');
            
            assert (isa (mbsys, 'mbdyn.pre.system'), ...
                'mbsys must be an mbdyn.pre.system object');
            
            if ~isempty (options.NEMOHSystem)
                assert (isa (options.NEMOHSystem, 'nemoh.simulation'), ...
                    'NEMOHSim must be a nemoh.simulation object' ); 
            end
            
            if isempty (options.MBDynInputFile)
                options.MBDynInputFile = fullfile (hsys.simu.caseDir, 'mbdyn_input_file.mbd');
            else
                assert (ischar (options.MBDynInputFile), ...
                    'MBDynInputFile must be string containing the file path where the MBDyn input file should be generated');
            end
            
            self.mBDynSystem = mbsys;
            self.hydroSystem = hsys;
            self.readyToRun = false;
            self.mBDynInputFile = options.MBDynInputFile;
            self.logger = wsim.logger ();
            
        end
        
        function addPTO (self, pto)
            
            assert (isa (pto, 'wsim.powertakeoff'), ...
                'pto must be a wsim.powertakeoff');
            
            self.readyToRun = false;
            
        end
        
        function prepare (self)
            
            self.mapPTOForceInds ();
            self.mapHydroForceInds ();
            
            % TODO: check if NEMOH data needs updated
            
            self.readyToRun = true;
            
        end
        
        function data = run (self, varargin)
            % Run the WEC simulation
            
            options.OutputFilePrefix = fullfile (self.hydroSystem.simu.caseDir, 'output', 'wsim');
            options.Verbosity = 0;
            
            options = parse_pv_pairs (options, varargin);
            
            check.isScalarInteger (options.Verbosity, true, 'Verbosity');
            
            assert (options.Verbosity >= 0, ...
                'Verbosity must an integer greater than or equal to zero');
            
            assert (ischar (options.OutputFilePrefix), ...
                'OutputFilePrefix must be a string');
            
            if self.readyToRun == false
                error ('Simulation is not ready to run, have you run ''prepare'' yet?');
            end
            
            % generate input file and start mbdyn
            
            % clear out previous files
            delete ([options.OutputFilePrefix, '.*']);
            
            % create the communicator object. As an mbdyn system object is
            % supplied, the mbdyn input file will be generated
            % automatically
            mb = mbdyn.mint.MBCNodal ('MBDynPreProc', mbsys, ...
                'UseMoments', true, ...
                'MBDynInputFile', self.mBDynInputFile, ...
                'OverwriteInputFile', true, ...
                'OutputPrefix', options.OutputFilePrefix ...
                );
            
            % ensure MBCNodal is destroyed in the event of a problem
            % (closes communication to MBDyn)
            CC = onCleanup (delete (mb));
            
            mb.start ('Verbosity', options.Verbosity);
            
            %%
            nnodes = mb.GetNodes ();
            
            time = mbsys.problems{1}.initialTime;
            ind = 1;
            
            status = mb.GetMotion ();
            
            if status ~= 0
                error ('mbdyn returned %d, aborting sim, check output file:\n%s\nfor clues at to why this happened.', status, mb.MBDynOutputFile)
            end
            
            eul = zeros (3,nnodes);
            
            R = mb.GetRot();
            for Rind = 1:size (R,3)
                om = mbdyn.pre.orientmat ('orientation', R(:,:,Rind));
                eul(1:3,Rind) = om.euler123();
            end
            
            pos = [ mb.NodePositions();
                    eul ];
            
            vel = [ mb.NodeVelocities();
                    mb.NodeOmegas() ];
            
            accel = [ mb.NodeAccelerations();
                      mb.NodeAngularAccels() ];
                  
            forces = zeros (3, nnodes);
            moments = zeros (3, nnodes);
            
            [hydroforces, out] = hsys.hydroForces (time, pos, vel, accel);
            
            for ind = 1:size (self.hydroNodeIndexMap, 1)
                forces (1:3,self.hydroNodeIndexMap(ind,1)) = ...
                    forces (1:3,self.hydroNodeIndexMap(ind,1)) + hydroforces(1:3,ind);
            end
            
            % PTO forces
            for ptoind = 1:size (self.ptoIndexMap, 1)
                
                ptoinfo = self.powerTakeOffs{ptoind}.loggingInfo ();
                
                tmpcell = cell(1,ptoinfo.NAvailable+1);
                
                [tmpcell{:}] =  self.powerTakeOffs{ptoind}.forceAndTorque ();
                
                ptoForceAndTorque = tmpcell{1};
                
                for outind = 1:ptoinfo.NAvailable
                    self.logger.logVal (ptoinfo.AvailableNames{outind}, 
                end
                
                FptoVec(1:3,1,ind) = F_pto(:,1);
                forces (1:3,1) = forces (1:3,1,ind) + F_pto(:,1);
                forces (1:3,2) = forces (1:3,2,ind) + F_pto(:,2);
                
            end

            
                
            % set the forces
            mb.F (forces(1:3,:));
            mb.M (forces(4:6,:));
            
            mbconv = mb.applyForcesAndMoments (false);
            
            data = initData
            
            % accept the last data into the time history of solutions
            hsys.advanceStep (time(1), vel(:,:,ind), accel(:,:,ind));
            
            ind = 2;
            
            plotvectors = false;
%             checkoutputs = false;
            miniters = 0;
            maxiters = self.mBDynSystem.problems{1}.maxIterations;
            absforcetol = 100;
            relforcetol = 1e-5;
            
            if plotvectors
                figure;
                hvectplotax = axes;
            end
            
            test_ptoforce(1) = 0;
            xRpto(1) = 0;
            vRpto(1) = 0;
            
            tic
            while status == 0
                
                status = mb.GetMotion ();
                
                if status ~= 0
                    continue;
                end
                
                time(ind) = time(ind-1) + self.mBDynSystem.problems{1}.timeStep;
                
                R = mb.GetRot();
                
                for Rind = 1:size (R,3)
                    om = mbdyn.pre.orientmat ('orientation', R(:,:,Rind));
                    eul(1:3,Rind) = om.euler123();
                end
                
                pos(:,:,ind) = [ mb.NodePositions(); eul];
                vel(:,:,ind) = [ mb.NodeVelocities(); mb.NodeOmegas() ];
                accel(:,:,ind) = [ mb.NodeAccelerations(); mb.NodeAngularAccels() ];
                
                [hydroforces, out] = hsys.hydroForces (time(ind), pos(:,:,ind), vel(:,:,ind), accel(:,:,ind));
                
                forces (:,:,ind) = hydroforces;
                
                % PTO force
                [F_pto, ptoforce(ind), test_xRpto(ind), test_vRpto(ind)] =  pto.force ();
                
                FptoVec(1:3,1,ind) = F_pto(:,1);
                forces (1:3,1,ind) = forces (1:3,1,ind) + F_pto(:,1);
                forces (1:3,2,ind) = forces (1:3,2,ind) + F_pto(:,2);
                
                mb.F (forces(1:3,:,ind));
                mb.M (forces(4:6,:,ind));
                
                mbconv = mb.applyForcesAndMoments (false);
                
                status = mb.GetMotion ();
                
                if status ~= 0
                    
                    F_ExcitLin(:,:,ind) = out.F_ExcitLin;
                    F_ViscousDamping(:,:,ind) = out.F_ViscousDamping;
                    F_AddedMass(:,:,ind) = out.F_AddedMass;
                    F_Restoring(:,:,ind) = out.F_Restoring;
                    F_RadiationDamping(:,:,ind) = out.F_RadiationDamping;
                    F_ExcitNonLin(:,:,ind) = out.F_ExcitNonLin;
                    F_MorrisonElement(:,:,ind) = out.F_MorrisonElement;
                    F_Excit(:,:,ind) = out.F_Excit;
                    F_ExcitRamp(:,:,ind) = out.F_ExcitRamp;
                    
                    ind = ind + 1;
                    
                    continue;
                end
                
                R = mb.GetRot();
                
                for Rind = 1:size (R,3)
                    om = mbdyn.pre.orientmat ('orientation', R(:,:,Rind));
                    eul(1:3,Rind) = om.euler123();
                end
                
                pos(:,:,ind) = [ mb.NodePositions(); eul];
                vel(:,:,ind) = [ mb.NodeVelocities(); mb.NodeOmegas() ];
                accel(:,:,ind) = [ mb.NodeAccelerations(); mb.NodeAngularAccels() ];
                
                % Hydrodynamic forces
                [newhydroforces, out] = hsys.hydroForces (time(ind), pos(:,:,ind), vel(:,:,ind), accel(:,:,ind));
                
                forces (:,:,ind) = newhydroforces;
                
                % PTO force
                [F_pto, test_ptoforce(ind), test_xRpto(ind), test_vRpto(ind)] =  pto.force ();
                
                FptoVec(1:3,1,ind) = F_pto(:,1);
                forces (1:3,1,ind) = forces (1:3,1,ind) + F_pto(:,1);
                forces (1:3,2,ind) = forces (1:3,2,ind) + F_pto(:,2);
                
                mb.F (forces(1:3,:,ind));
                mb.M (forces(4:6,:,ind));
                
                mbconv = mb.applyForcesAndMoments (false);
                
                
                forcediff = abs (hydroforces - newhydroforces);
                
                maxforces = max(hydroforces, newhydroforces);
                relforcediff = abs(forcediff) ./ abs(maxforces);
                relforcediff(maxforces == 0) = 0;
                itercount = 1;
                while mbconv ~= 0 ...
                        || itercount < miniters ...
                        || (max (forcediff(:)) > absforcetol) ...
                        || (ind > 3 && (max (relforcediff(:)) > relforcetol))
                    
                    % store the previously calculated hydrodynamic forces
                    hydroforces = newhydroforces;
                    
                    status = mb.GetMotion ();
                    
                    if status ~= 0
                        break;
                    end
                    
                    R = mb.GetRot();
                    
                    for Rind = 1:size (R,3)
                        om = mbdyn.pre.orientmat ('orientation', R(:,:,Rind));
                        eul(1:3,Rind) = om.euler123();
                    end
                    
                    pos(:,:,ind) = [ mb.NodePositions(); eul];
                    vel(:,:,ind) = [ mb.NodeVelocities(); mb.NodeOmegas() ];
                    accel(:,:,ind) = [ mb.NodeAccelerations(); mb.NodeAngularAccels() ];
                    
                    % Hydrodynamic forces
                    [newhydroforces, out] = hsys.hydroForces (time(ind), pos(:,:,ind), vel(:,:,ind), accel(:,:,ind));
                    
                    forces (:,:,ind) = newhydroforces;
                    
                    % PTO force
                    [F_pto, ptoforce(ind), test_xRpto(ind) vRpto(ind)] = pto.force ();
                    
                    FptoVec(1:3,1,ind) = F_pto(:,1);
                    forces (1:3,1,ind) = forces (1:3,1,ind) + F_pto(:,1);
                    forces (1:3,2,ind) = forces (1:3,2,ind) + F_pto(:,2);
                    
                    mb.F (forces(1:3,:,ind));
                    mb.M (forces(4:6,:,ind));
                    
                    mbconv = mb.applyForcesAndMoments (false);
                    
                    itercount = itercount + 1;
                    
                    if itercount > maxiters
                        error ('mbdyn iterations exceeded max allowed');
                    end
                    
                    forcediff = abs (hydroforces - newhydroforces);
                    maxforces = max(hydroforces, newhydroforces);
                    relforcediff = abs(forcediff) ./ abs(maxforces);
                    relforcediff(maxforces == 0) = 0;
                    
                end
                
                status = mb.GetMotion ();
                
                if status ~= 0
                    
                    F_ExcitLin(:,:,ind) = out.F_ExcitLin;
                    F_ViscousDamping(:,:,ind) = out.F_ViscousDamping;
                    F_AddedMass(:,:,ind) = out.F_AddedMass;
                    F_Restoring(:,:,ind) = out.F_Restoring;
                    F_RadiationDamping(:,:,ind) = out.F_RadiationDamping;
                    F_ExcitNonLin(:,:,ind) = out.F_ExcitNonLin;
                    F_MorrisonElement(:,:,ind) = out.F_MorrisonElement;
                    F_Excit(:,:,ind) = out.F_Excit;
                    F_ExcitRamp(:,:,ind) = out.F_ExcitRamp;
                    
                    ind = ind + 1;
                    
                    break;
                end
                
                mb.F (forces(1:3,:,ind));
                mb.M (forces(4:6,:,ind));
                
                mbconv = mb.applyForcesAndMoments (true);
                
                F_ExcitLin(:,:,ind) = out.F_ExcitLin;
                F_ViscousDamping(:,:,ind) = out.F_ViscousDamping;
                F_AddedMass(:,:,ind) = out.F_AddedMass;
                F_Restoring(:,:,ind) = out.F_Restoring;
                F_RadiationDamping(:,:,ind) = out.F_RadiationDamping;
                F_ExcitNonLin(:,:,ind) = out.F_ExcitNonLin;
                F_MorrisonElement(:,:,ind) = out.F_MorrisonElement;
                F_Excit(:,:,ind) = out.F_Excit;
                F_ExcitRamp(:,:,ind) = out.F_ExcitRamp;
                
                % accept the last data into the time history of solutions
                hsys.advanceStep (time(ind), vel(:,:,ind), accel(:,:,ind));
                
                ind = ind + 1;
                
            end
            
            [F_Total, F_AddedMassCorrected] = correctAddedMassForce (hsys, forces, F_AddedMass, accel);
            toc;
            
            fprintf (1, 'Reached time %f, in %d steps\n', time(end), ind-1);
            
        end
        
    end
    
    methods (Access = private)
        
        function [forces, moments] = applyPTOForces (self, forces, moments)
            
                       % PTO forces
            for ind = 1:size (self.ptoIndexMap, 1)
                
                if isa (pto, 'mbdyn.mint.twoNodeTranslationalForce')
                    
                    [F_pto, ptoforce(ind), test_xRpto(ind), test_vRpto(ind)] =  pto.force ();

                    forces (1:3,self.ptoIndexMap(ind,1)) = ...
                        forces (1:3,self.ptoIndexMap(ind,1)) + F_pto(1:3,1);

                    forces (1:3,self.ptoIndexMap(ind,1)) = ...
                        forces (1:3,self.ptoIndexMap(ind,2)) + F_pto(1:3,2);
                
                elseif isa (pto, 'mbdyn.mint.twoNodeRotationalForce')
                    
                    [Torque_pto, ptoforce(ind), test_xRpto(ind), test_vRpto(ind)] =  pto.force ();

                    forces (1:3,self.ptoIndexMap(ind,1)) = ...
                        forces (1:3,self.ptoIndexMap(ind,1)) + Torque_pto(1:3,1);

                    forces (1:3,self.ptoIndexMap(ind,1)) = ...
                        forces (1:3,self.ptoIndexMap(ind,2)) + Torque_pto(1:3,2);
                    
                else
                    error ('Unsupported pto force type');
                end
                
            end
            
            
        end
        
        function mapPTOForceInds (self)
            
            strinfo = self.mBDynSystem.externalStructuralInfo ();
            
            self.ptoIndexMap = zeros (numel (self.powerTakeOffs), 2);
            
            for ptoind = 1:numel(self.powerTakeOffs)
                
                foundptonode = false;
                for nodeind = 1:numel (strinfo.Nodes)
                    
                    if strinfo.Nodes{nodeind} == self.powerTakeOffs{ptoind}.referenceNode
                        
                        foundptonode = true;
                        
                        self.ptoIndexMap(ptoind, 1) = nodeind;
                        
                        break;
                    
                    end
                
                end
                
                if foundptonode == false
                    error ('Could not find reference node for PTO %d in MBDyn system', ptoind);
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
                    error ('Could not find non-reference node for PTO %d in MBDyn system', ptoind);
                end
                
            end
            
        end
        
        function mapHydroForceInds (self)
            
            strinfo = self.mBDynSystem.externalStructuralInfo ();
            
            self.hydroNodeIndexMap = zeros (numel (self.hydroSystem.bodyMBDynNodes), 2);
            
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
                    error ('Could not find corresponding external structural node for hydro node %d', hydronodeind);
                end
                
            end
            
        end
        

        function initDataStructures ( self, storepos, ...
                                            storevel, ... 
                                            storeaccel, ...
                                            storehydroforce, ...
                                            storeexcitforce, ...
                                            storeexcitlinforce, ...
                                            storeexcitnonlinforce, ...
                                            storeradforce, ...
                                            storerestoringforce, ...
                                            storevisdampingforce, ...
                                            storeptoforce, ...
                                            storelinptoxRpto )
                        
            
            extforceinfo = self.mBDynSystem.externalStructuralInfo ();

            % preallocate data vectors
            nsteps = (simu_wsim.endTime - simu_wsim.startTime) ./ simu_wsim.dt + 1;
            
            % always store time, as is is independent variabl for other
            % vars in logger object
            self.logger.addVariable ( 'Time', [1, 1], ...
                                      'Desc', 'main hydrodynamic/multibody time step', ...
                                      'Pre', nsteps );
            
            if storepos
                self.logger.addVariable ( 'Positions', [6, extforceinfo.NNodes], ...
                                          'Desc', 'cartesian and angular positions of all structural external nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if storevel
                self.logger.addVariable ( 'Velocities', [6, extforceinfo.NNodes], ...
                                          'Desc', 'cartesian and angular velocities of all structural external nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if storeaccel
                self.logger.addVariable ( 'Accelerations', [6, extforceinfo.NNodes], ...
                                          'Desc', 'cartesian and angular accelerations of all structural external nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if storehydroforce
                self.logger.addVariable ( 'ForceHydro', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'summ of all hydrodynamic forces for all hydro nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if storeexcitforce
                self.logger.addVariable ( 'ForceExcitation', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'sum of linear and nonlinear hydrodynamic excitation force for all hydro nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if storeexcitlinforce
                self.logger.addVariable ( 'ForceExcitationLin', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'linear hydrodynamic excitation forces for all hydro nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if storeexcitnonlinforce
                self.logger.addVariable ( 'ForceExcitationNonLin', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'nonlinear hydrodynamic excitation forces for all hydro nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if storeradforce
                self.logger.addVariable ( 'ForceRadiationDamping', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'hydrodynamic radiation and damping forces for all hydro nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if storerestoringforce
                self.logger.addVariable ( 'ForceRestoring', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'hydrodynamic restoring forces for all hydro nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if storemorrisonforce
                self.logger.addVariable ( 'ForceMorrison', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'morrison forces for all hydro nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if storevisdampingforce
                self.logger.addVariable ( 'ForceViscousDmping', [6, self.hydroSystem.nHydroBodies], ...
                                          'Desc', 'hydrodynamic viscous damping forces for all hydro nodes', ...
                                          'Pre', nsteps, ...
                                          'Indep', 'Time' );
            end
            
            if storeptoforce
                data.PTOForces = cell (size (self.powerTakeOffs));
                for ptoind = 1:numel (data.PTOForces)
                    n = self.powerTakeOffs{ind}.forceSize ();
                    data.PTOForces{ptoind} = nan ([n, 1, nsteps]);
                end
            end
            
            if storelinptoxRpto
                data.Positions = nan ([1, 1, nsteps]);
            end
            
            xRpto = ptoforce;
            vRpto = ptoforce;
            xRptoVec = repmat ([0;0;0], [1, nsteps]);
            vRptoVec = repmat ([0;0;0], [1, nsteps]);
            FptoVec = repmat ([0;0;0], [1, 1, nsteps]);
            % twoNodeTranslationalForce test variables
            test_xRpto = xRpto;
            test_vRpto = vRpto;
            test_FptoVec = FptoVec;
            test_xRptoVec = xRptoVec;
            test_vRptoVec = vRptoVec;
            
            
            F_ViscousDamping = repmat (out.F_ViscousDamping, [1,1,nsteps]);
            F_AddedMass = repmat (out.F_AddedMass, [1,1,nsteps]);
            F_Restoring = repmat (out.F_Restoring, [1,1,nsteps]);
            F_RadiationDamping = repmat (out.F_RadiationDamping, [1,1,nsteps]);
            F_ExcitNonLin = repmat (out.F_ExcitNonLin, [1,1,nsteps]);
            F_MorrisonElement = repmat (out.F_MorrisonElement, [1,1,nsteps]);
            F_Excit = repmat (out.F_Excit, [1,1,nsteps]);
            F_ExcitRamp = repmat (out.F_ExcitRamp, [1,1,nsteps]);
            
            
        end
        
    end
    
end