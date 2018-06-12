classdef hydroSystem < handle
% class representing a collection of hydrodynamically interacting bodies
%
% Description
%
% wsim.hydroSystem is a class used to simulate a one or more bodies
% interaction with fluid waves and also, optionally, body-to-body
% interaction.
%
% Methods:
%
%  addHydroBodies
%  initialiseHydrobodies
%
% 
% See also: wsim.hydroBody, wsim.wecSim
%
    
    properties (GetAccess = public, SetAccess = private)
        
        nHydroBodies;
        odeSimInitialised = false;
        hydroBodiesInitialised = false;
        simu;        % simulation setings (wsim.simSettings object)
        waves;       % wave settings (wsim.waveSettings object)
        bodyMBDynNodes;
        caseDirectory;
        
    end
    
    properties (GetAccess = private, SetAccess = private)
        
        
        hydroBodies; % array of wsim.hydroBody objects
        hydroBodyInds = []; % array of indices of the hydrobodies
        
    end
    
    
    methods
        
        function self = hydroSystem (waves, simu, hydrobodies)
            % wsim.hydroSystem class constructor
            %
            % Syntax
            %
            % hsys = hydroSystem (waves, simu)
            % hsys = hydroSystem (..., hydrobodyfiles)
            % 
            % Description
            %
            % wsim.hydroSystem is a class used to simulate a one or more
            % bodies interaction with fluid waves and also, optionally,
            % body-to-body interaction.
            %
            % Input
            %
            %  waves - WEC-Sim waves class
            %
            %  simu - WEC-Sim simu class
            %
            %  hydrobodies - optional array of hydrobody objects.
            %    Additional bodies can be added later with the
            %    'addHydroBodies' method.
            %
            % Output
            %
            %  hsys - hydroSystem object
            % 
            
            
            self.simu = simu;
            self.waves = waves;
            
            self.odeSimInitialised = false;
            
            % Waves and Simu: check inputs
            self.waves.checkInputs ();
            
            % initialise the number of WEC bodies in the sim to 0, this
            % will be incremented as they are added
            self.simu.numWecBodies = 0;
            
            % get the case directory from the sim settings
            self.caseDirectory = self.simu.caseDir;
            
            % add the supplied hyrobodies to the hydrodynamic system
            addHydroBodies (self, hydrobodies);
            
        end
        
        function addHydroBodies (self, hydrobodies)
            % add one or more hydrodynamic bodies to the system
            %
            % Syntax
            %
            %  addHydroBodies (hsys, hydrobodyfiles)
            %
            %
            % Input
            %
            %  hsys - wsim.hydroSystem object
            %
            %  hydrobodies - array of one or more hydrobody objects to be
            %    added to the system
            %
            
            
            if self.odeSimInitialised
                
                error ('An ode simulation has already been initialised so you cannot add new bodies. To add new bodies first call the timeDomainSimReset method.')
                
            else
                
                if ~isa (hydrobodies, 'wsim.hydroBody')
                    error ('hydrobodyfiles must be a hydrobody object, or array of hydrobody objects');
                end

                if isempty (self.hydroBodyInds)
                    nextbodyind = 1;
                else
                    nextbodyind = self.hydroBodyInds(end) + 1;
                end

                for ind = 1:numel (hydrobodies)

                    % update the list of body indices
                    self.hydroBodyInds = [ self.hydroBodyInds, nextbodyind ];
                    
                    hydrobodies(ind).bodyNumber = nextbodyind;
                    
                    % the following complication is necessary to support
                    % Octave, otherwise we could just do
                    % self.hydroBodies = [self.hydroBodies, hydrobodies(ind)];
                    if isempty (self.hydroBodies)
                        self.hydroBodies = hydrobodies(ind);
                    else
                        % add the new hydrobodies to the collection
                        self.hydroBodies(end+1) = hydrobodies(ind);
                    end
                
                    % increment the number of bodies in simu
                    self.simu.numWecBodies = self.simu.numWecBodies + 1;
                    
                    nextbodyind = nextbodyind + 1;

                end
                
                self.hydroBodiesInitialised = false;
            
            end
            
        end
        
        function initialiseHydrobodies (self, varargin)
            % performs some organisation tasks for the system
            %
            % Syntax
            %
            % initialiseHydrobodies (hsys)
            % initialiseHydrobodies (..., 'Parameter', Value)
            %
            % Description
            %
            % initialiseHydrobodies performs some preprocessing task to
            % organise the system which are necessary before proceeding to
            % set up a tme domain system. It should be called before
            % calling the timeDomainSimSetup method.
            %
            % Input
            %
            %  hsys - wsim.hydroSystem object
            %
            % Addtional arguments may be supplied as parameter-value pairs.
            % The available options are:
            %
            %  'MultiConditionRun' - 
            %
            %  'MultiConditionRunIndex' - 
            %
            % Output
            %
            %
            %
            % See Also: wsim.hydroSystem.timeDomainSimSetup
            %
            
            options.MultiConditionRun = false;
            options.MultiConditionRunIndex = 1;
            
            options = parse_pv_pairs (options, varargin);
            
            % Bodies: count, check inputs, read hdf5 file
            
%             for bodyind = 1:numel (self.hydroBodies)
% 
%                 if self.hydroBodies(bodyind).nhBody == 0
% 
%                     numHydroBodies = numHydroBodies + 1;
% 
%                 else
% 
%                     self.hydroBodies(bodyind).bemioFlag = 0;
%                     self.hydroBodies(bodyind).massCalcMethod = 'user';
% 
%                 end
%             end

            % Read Inputs for Multiple Conditions Run 
            if options.MultiConditionRun == true
                
                for n = 1:length (mcr.cases(1,:))
                    
                    if iscell (mcr.cases)
                        eval ([mcr.header{n} '= mcr.cases{imcr,n};']);
                    else
                        eval ([mcr.header{n} '= mcr.cases(imcr,n);']);
                    end
                end
                
            end
            
            for bodyind = 1:numel (self.hydroBodies)
                
%                 self.hydroBodies(bodyind).initialiseSimParams (waves, simu);

                self.hydroBodies(bodyind).setCaseDirectory (self.caseDirectory);
                
                self.hydroBodies(bodyind).checkInputs ();
                
                % Determine if hydro data needs to be reloaded from h5
                % file, or if hydroData was stored in memory from a
                % previous run.
                if options.MultiConditionRun == true ...
                        && self.simu.reloadH5Data == 0 ...
                        && options.MultiConditionRunIndex > 1 
                    
                    self.hydroBodies(bodyind).loadHydroData (hydroData(bodyind));
                    
                else
                    % load the hydrodynamic data from file
                    self.hydroBodies(bodyind).loadHydroData ();
                    
%                     % check water depth is the same as the previously added
%                     % body
%                     if bodyind > 1
%                         if self.hydroBodies(bodyind).hydroData.simulation_parameters.water_depth ...
%                                 ~= self.hydroBodies(bodyind-1).hydroData.simulation_parameters.water_depth
%                             
%                             error ('Water depth for body %d is not the same as water depth for body %d', bodyind, bodyind-1)
%                             
%                         end
%                     end
                    
                end
                
                self.hydroBodies(bodyind).bodyTotal = self.simu.numWecBodies;
                
                if self.simu.b2b == 1
                    self.hydroBodies(bodyind).lenJ = zeros (6*self.hydroBodies(bodyind).bodyTotal, 1);
                else
                    self.hydroBodies(bodyind).lenJ = zeros (6, 1);
                end
                
            end
            
            self.hydroBodiesInitialised = true;
            
        end
        
        function timeDomainSimSetup (self)
            % prepare wsim.hydroSystem for a time domain simulation
            %
            % Syntax
            %
            %  timeDomainSimSetup (hsys)
            %
            % Description
            %
            % timeDomainSimSetup prepares the hydrodynamic system to
            % perform a time domain simulation of the system.
            %
            % Input
            %
            %  hsys - wsim.hydroSystem object
            %
            
            % simulation setup
            self.simu.checkInputs ();
            self.simu.setupSim ();
            
            % wave setup
            self.waves.waveSetup( self.hydroBodies(1).hydroData.simulation_parameters.w, ...
                                  self.hydroBodies(1).hydroData.simulation_parameters.water_depth, ...
                                  self.simu.rampT, ...
                                  self.simu.dt, ...
                                  self.simu.maxIt, ...
                                  self.simu.g, ...
                                  self.simu.endTime );
                              
            % Non-linear hydro
            if (self.simu.nlHydro > 0) || (self.simu.paraview == 1)
                for bodyind = 1:length(self.hydroBodies(1,:))
                    self.hydroBodies(bodyind).bodyGeo (body(bodyind).geometryFile)
                end
            end
            
            % call each body's timeDomainSimSetup method
            for bodyind = 1:numel(self.hydroBodies)
                
                self.hydroBodies(bodyind).timeDomainSimSetup (self.waves, self.simu, self.hydroBodyInds(bodyind));
                
            end
            
            for bodyind = 1:numel(self.hydroBodies)
                self.hydroBodies(bodyind).adjustMassMatrix ();
            end
            
            self.odeSimInitialised = true;
            
        end
        
        function [mbnodes, mbbodies, mbelements] = makeMBDynComponents (self)
            % generates mbdyn preprocessor objects for the system
            %
            % Syntax
            %
            % [mbnodes, mbbodies, mbelements] = makeMBDynComponents (hsys)
            %
            % Description
            %
            % Generates mbdyn preprocessor objects for the system
            %
            % Input
            %
            %  hsys - wsim.hydroSystem object
            %
            % Output
            %
            %  mbnodes - cell array of mbdyn.pre.structuralNode6dof
            %   objects, one for each body in the system
            %
            %  mbbodies - cell array of mbdyn.pre.body objects, one for
            %   each body in the system
            %
            %  mbelements - cell array of other mbdyn preprocessor objects
            %   required for the system. The exact contents of this depends
            %   on the simulation settings
            %
            %
            
            if self.hydroBodiesInitialised
                            
                mbnodes = {};
                mbbodies = {};
                mbelements = {};
                
                input_list = {};
                
                % make the structural nodes and bodies
                for bodyind = 1:numel (self.hydroBodies)
                    
                    [node, body] = self.hydroBodies(bodyind).makeMBDynComponents ();
                
                    mbnodes = [mbnodes, {node}];
                    mbbodies = [mbbodies, {body}];
                    
                end
                
                self.bodyMBDynNodes = mbnodes;
                
                absnodes = {};
                forces = {};
                if self.simu.ssCalc == 2
                    % need to create abstract nodes for the moments of the
                    % structural nodes as these cannot be addressed
                    % directly
                    
                    for mbnodeind = 1:numel(mbnodes)
                        if isa (mbnodes{mbnodeind}, 'mbdyn.pre.structuralNode6dof')
                            % create three abstract nodes, one for each
                            % moment
                            
                            newabsnodes = {mbdyn.pre.abstractNode('AlgebraicOrDifferential', 'algebraic'), ...
                                           mbdyn.pre.abstractNode('AlgebraicOrDifferential', 'algebraic'), ...
                                           mbdyn.pre.abstractNode('AlgebraicOrDifferential', 'algebraic'), ...
                                           mbdyn.pre.abstractNode('AlgebraicOrDifferential', 'algebraic'), ...
                                           mbdyn.pre.abstractNode('AlgebraicOrDifferential', 'algebraic'), ...
                                           mbdyn.pre.abstractNode('AlgebraicOrDifferential', 'algebraic') };
                                       
                            absnodes = [ absnodes, newabsnodes ];
                            
                            % next create a moment which uses the output of
                            % these nodes to apply a couple to the
                            % structural node
                            drivecallers = { mbdyn.pre.nodeDrive(newabsnodes{1}, mbdyn.pre.linearDrive (0, -1)), ...
                                        	 mbdyn.pre.nodeDrive(newabsnodes{2}, mbdyn.pre.linearDrive (0, -1)), ...
                                             mbdyn.pre.nodeDrive(newabsnodes{3}, mbdyn.pre.linearDrive (0, -1)), ...
                                             mbdyn.pre.nodeDrive(newabsnodes{4}, mbdyn.pre.linearDrive (0, -1)), ...
                                        	 mbdyn.pre.nodeDrive(newabsnodes{5}, mbdyn.pre.linearDrive (0, -1)), ...
                                             mbdyn.pre.nodeDrive(newabsnodes{6}, mbdyn.pre.linearDrive (0, -1)) };
                            
                            fdc = mbdyn.pre.componentTplDriveCaller (drivecallers(1:3));
                            mdc = mbdyn.pre.componentTplDriveCaller (drivecallers(4:6));
                            
                            forces = [forces, { mbdyn.pre.structuralForce(mbnodes{mbnodeind}, 'absolute', 'null', fdc), ...
                                                mbdyn.pre.structuralCouple(mbnodes{mbnodeind}, 'absolute', mdc) } ];
                            
                        end
                    end
                end
                
                for bodyind = 1:numel (self.hydroBodies)
                    
                    if self.simu.ssCalc == 2
                
                        if self.simu.b2b
                            % inputs are all the nodes velocities
                            for mbnodeind = 1:numel(mbnodes)
                                % forces on node
                                input_list = [input_list, {mbdyn.pre.nodeDrive(mbnodes{mbnodeind}, mbdyn.pre.directDrive (), 'String', 'XP[1]')}];
                                input_list = [input_list, {mbdyn.pre.nodeDrive(mbnodes{mbnodeind}, mbdyn.pre.directDrive (), 'String', 'XP[2]')}];
                                input_list = [input_list, {mbdyn.pre.nodeDrive(mbnodes{mbnodeind}, mbdyn.pre.directDrive (), 'String', 'XP[3]')}];
                                % moments on node, cannot be addressed
                                % directly, so use Node Drive
                                input_list = [input_list, {mbdyn.pre.nodeDrive(mbnodes{mbnodeind}, mbdyn.pre.directDrive (), 'String', 'Omega[1]')}];
                                input_list = [input_list, {mbdyn.pre.nodeDrive(mbnodes{mbnodeind}, mbdyn.pre.directDrive (), 'String', 'Omega[2]')}];
                                input_list = [input_list, {mbdyn.pre.nodeDrive(mbnodes{mbnodeind}, mbdyn.pre.directDrive (), 'String', 'Omega[3]')}];
                            end
                        else
                            % inputs are just this body's node velocities
                            input_list = { mbdyn.pre.nodeDrive(mbnodes{bodyind}, mbdyn.pre.directDrive (), 'String', 'XP[1]'), ...
                                           mbdyn.pre.nodeDrive(mbnodes{bodyind}, mbdyn.pre.directDrive (), 'String', 'XP[2]'), ...
                                           mbdyn.pre.nodeDrive(mbnodes{bodyind}, mbdyn.pre.directDrive (), 'String', 'XP[3]'), ...
                                           mbdyn.pre.nodeDrive(mbnodes{bodyind}, mbdyn.pre.directDrive (), 'String', 'Omega[1]'), ...
                                           mbdyn.pre.nodeDrive(mbnodes{bodyind}, mbdyn.pre.directDrive (), 'String', 'Omega[2]'), ...
                                           mbdyn.pre.nodeDrive(mbnodes{bodyind}, mbdyn.pre.directDrive (), 'String', 'Omega[3]') };
                        end

%                         output_node_list = mbdyn.pre.nodeDOF (mbnodes{bodyind}, 'DOFNumber', 1);
%                         output_node_list(2) = mbdyn.pre.nodeDOF (mbnodes{bodyind}, 'DOFNumber', 2);
%                         output_node_list(3) = mbdyn.pre.nodeDOF (mbnodes{bodyind}, 'DOFNumber', 3);
                        output_node_list = mbdyn.pre.nodeDOF (absnodes{(bodyind-1)*6+1}, 'AlgebraicOrDifferential', 'algebraic');
                        output_node_list(2) = mbdyn.pre.nodeDOF (absnodes{(bodyind-1)*6+2}, 'AlgebraicOrDifferential', 'algebraic');
                        output_node_list(3) = mbdyn.pre.nodeDOF (absnodes{(bodyind-1)*6+3}, 'AlgebraicOrDifferential', 'algebraic');
                        output_node_list(4) = mbdyn.pre.nodeDOF (absnodes{(bodyind-1)*6+4}, 'AlgebraicOrDifferential', 'algebraic');
                        output_node_list(5) = mbdyn.pre.nodeDOF (absnodes{(bodyind-1)*6+5}, 'AlgebraicOrDifferential', 'algebraic');
                        output_node_list(6) = mbdyn.pre.nodeDOF (absnodes{(bodyind-1)*6+6}, 'AlgebraicOrDifferential', 'algebraic');

%                         sys = ss( self.hydroBodies(bodyind).hydroForce.ssRadf.A, ...
%                                   self.hydroBodies(bodyind).hydroForce.ssRadf.B, ...
%                                   self.hydroBodies(bodyind).hydroForce.ssRadf.C, ...
%                                   self.hydroBodies(bodyind).hydroForce.ssRadf.D );
%                               
%                         sys = prescale (sys, { self.hydroBodies(bodyind).hydroData.simulation_parameters.w(1), ...
%                                                self.hydroBodies(bodyind).hydroData.simulation_parameters.w(end) });
%               
%                         if all (all (sys.D == 0))
%                             mbelements = [ mbelements, ...
%                                            { mbdyn.pre.stateSpaceMIMO( size (self.hydroBodies(bodyind).hydroForce.ssRadf.A, 1), ...
%                                                                 sys.A, ...
%                                                                 sys.B, ...
%                                                                 sys.C, ...
%                                                                 output_node_list, ...
%                                                                 input_list, ...
%                                                                 'Balance', 'yes') }, ...
%                                             ];
%                         else
%                             mbelements = [ mbelements, ...
%                                            { mbdyn.pre.stateSpaceMIMO( size (self.hydroBodies(bodyind).hydroForce.ssRadf.A, 1), ...
%                                                                 sys.A, ...
%                                                                 sys.B, ...
%                                                                 sys.C, ...
%                                                                 output_node_list, ...
%                                                                 input_list, ...
%                                                                 'Balance', 'yes', ...
%                                                                 'D', sys.D ) }, ...
%                                             ];
%                         end

                        if all (all (sys.D == 0))
                            mbelements = [ mbelements, ...
                                           { mbdyn.pre.stateSpaceMIMO( size (self.hydroBodies(bodyind).hydroForce.ssRadf.A, 1), ...
                                                                self.hydroBodies(bodyind).hydroForce.ssRadf.A, ...
                                                                self.hydroBodies(bodyind).hydroForce.ssRadf.B, ...
                                                                self.hydroBodies(bodyind).hydroForce.ssRadf.C, ...
                                                                output_node_list, ...
                                                                input_list, ...
                                                                'Balance', 'yes') }, ...
                                            ];
                        else
                            mbelements = [ mbelements, ...
                                           { mbdyn.pre.stateSpaceMIMO( size (self.hydroBodies(bodyind).hydroForce.ssRadf.A, 1), ...
                                                                self.hydroBodies(bodyind).hydroForce.ssRadf.A, ...
                                                                self.hydroBodies(bodyind).hydroForce.ssRadf.B, ...
                                                                self.hydroBodies(bodyind).hydroForce.ssRadf.C, ...
                                                                output_node_list, ...
                                                                input_list, ...
                                                                'Balance', 'yes', ...
                                                                'D', self.hydroBodies(bodyind).hydroForce.ssRadf.D ) }, ...
                                            ];
                        end
                                    
                    end
                end
                
                mbelements = [ mbelements, forces ];
                mbnodes = [mbnodes, absnodes];
                
            else
                error ('You must call initialiseHydrobodies before attempting to create the mbdyn components.');
            end
            
            
        end
        
        function [forces, breakdown] = hydroForces (self, t, pos, vel, accel)
            % calculates forces acting on all bodies in a hydrodynamic system
            %
            % [forces, breakdown] = hydroForces (hsys, t, x, vel, accel)
            %
            % Description
            %
            % hydroForces calculates the hydrodynamic forces acting on all
            % bodies in a hydrodynamic system.
            %
            % Input
            %
            %  hsys - hydroSystem object
            % 
            %  t - current time
            %
            %  pos - (6 x nBodies) cartesian and angular positions of the
            %    bodies in the system. Each column of pos is the positions
            %    of one body, with the first three elements being the x, y,
            %    and z positions and the second the three the angular
            %    orientations in radians using the extrinsic Euler 123
            %    (xyz) convention.
            %
            %  vel - (6 x nBodies) cartesian and angular velocities of the
            %    bodies in the system. The velocites correspond to the
            %    equvalent values in 'pos'.
            %
            %  accel - (6 x nBodies) cartesian and angular accelerations of
            %    the bodies in the system. The accelerations correspond to
            %    the equvalent values in 'pos'.
            %
            % Output
            %
            %  forces - (6 x nBodies) matrix of translational and angular
            %    forces acting on the body. 
            %
            %  breakdown - structure containing a more detailed breakdown
            %    of the forces and additional information. The fields which
            %    will be present depend on the details of the simulation
            %    chosen, but can include the following:
            %
            %    F_ExcitLin : (6 x nBodies) matrix of linear wave
            %      excitation forces. See also F_Excit and F_ExcitRamp
            %      below.
            %
            %    F_ExcitNonLin : (6 x nBodies) matrix of non-linear wave
            %      excitation forces (if calculated). Will be all zeros if
            %      non-linear excitation forces have not been activated for
            %      the simulation. See also F_Excit and F_ExcitRamp below.
            %
            %    F_Excit : (6 x nBodies) matrix of total wave
            %      excitation forces. This is the sum of F_ExcitLin and
            %      F_ExcitNonLin. If a ramp function is applied, the actual
            %      forces applied to the bodies are in F_ExcitRamp
            %
            %    F_ExcitRamp : (6 x nBodies) matrix of wave excitation
            %      forces (i.e. as in F_Excit), but with a ramp function
            %      applied. These are the actual excitation forces applied
            %      to the bodies during the simulation.
            %  
            %    F_ViscousDamping : (6 x nBodies) matrix of viscous damping
            %      forces.
            %
            %    F_AddedMass : (6 x nBodies) matrix of added mass forces.
            %      Note that after simulation these should be modified
            %      using the correctAddedMassForce method, to get the
            %      'real' added mass forces.
            %
            %    F_RadiationDamping : (6 x nBodies) matrix of wave
            %      radiation damping forces.
            %
            %    F_Restoring : (6 x nBodies) matrix of restoring (buoyancy)
            %      forces.
            %
            %    F_MorrisonElement : (6 x nBodies) matrix of morrison
            %      element forces.
            %
            %    BodyHSPressure : 
            %
            %    WaveNonLinearPressure : The Froude–Krylov pressure
            %
            %    WaveLinearPressure : 
            %
            
            forces = nan * ones (6, self.simu.numWecBodies );
            
            for bodyind = 1:numel(self.hydroBodies)
                [forces(:,bodyind), bodybreakdown] ...
                        = hydroForces ( self.hydroBodies(bodyind), ...
                                        t, ...
                                        pos(:,bodyind), ...
                                        vel, ...
                                        accel );
                
                if nargout > 1
                    if bodyind == 1
                        fnames = fieldnames (bodybreakdown);
                    end

                    for bdfind = 1:numel(fnames)
                        if isempty (bodybreakdown.(fnames{bdfind}))
                            breakdown.(fnames{bdfind}) = [];
                        else
                            breakdown.(fnames{bdfind})(:,bodyind) = bodybreakdown.(fnames{bdfind});
                        end
                    end
                end
            end
            
        end
        
        function advanceStep (self, t, vel, accel)
            % advance to the next time step, store data as required
            %
            % Syntax
            %
            % advanceStep (hsys, t, vel, accel)
            %
            % Description
            %
            % advanceStep advances the solution to the next time step,
            % accepting the current time step solution and data into stored
            % solution histories. 
            %
            % Input
            %
            %  hsys - wsim.hydroSystem object
            %
            %  t - the current simulation time
            %
            %  vel - (6 x nbodies) matrix of body velocities and angular
            %    velocities at time t.
            %
            %  accel - (6 x nbodies) matrix of body accelerations and
            %    angular accelerations at time t.
            %
            % Output
            %
            %  none
            %
            %

            for bodyind = 1:numel(self.hydroBodies)
                if self.simu.b2b
                    self.hydroBodies(bodyind).advanceStep (t, vel(:), accel(:));
                else
                    self.hydroBodies(bodyind).advanceStep (t, vel(:,bodyind), accel(:,bodyind));
                end
            end

        end
        
        function [F_Total, F_AddedMass] = correctAddedMassForce (self, forceTotal, forceAddedMass, accel)
            % recalculate the added mass and total forces on bodies
            %
            % Syntax
            %
            % [F_Total, F_AddedMass] = correctAddedMassForce (hsys, forceTotal, forceAddedMass, accel)
            %
            % Description
            %
            % 
            %
            % Input
            %
            %  hsys - hydroSystem object
            %
            %  forceTotal - the total force calculated in a transient
            %    simulation before added mass force correction (as output
            %    by hydroSystem.hydroForces)
            %
            %  forceAddedMass - the added mass force calculated in a
            %    transient simulation before added mass force correction
            %    (as output by hydroSystem.hydroForces)
            %
            %  accel - (6 x n x nBodies) matrix of accelerations calculated
            %    during a transient simulation 
            %
            % Output
            %
            %  F_Total - the corrected total forces for each node (taking
            %   account the corected 
            %
            %  F_AddedMass - the corected added mass forces
            %
            %
            
            F_Total = forceTotal + forceAddedMass;
            F_AddedMass = nan * ones (size (forceAddedMass));
            
            for bodyind = 1:numel(self.hydroBodies)
                self.hydroBodies(bodyind).restoreMassMatrix ();
                
                % body.forceAddedMass returns an (n x 6) matrix of values,
                % with the rows being the time series history. We reshape
                % this into a (6 x 1 x n) matrix for insertion into the
                % F_AddedMass matrix in the appropriate location. To get
                % the right shape we have to transpose the result of
                % body.forceAddedMass before reshaping
                F_AddedMass(:,bodyind,:) = reshape ( ...
                    self.hydroBodies(bodyind).forceAddedMass(squeeze(accel(:,bodyind,:)).', self.simu.b2b).', ...
                                                     6, 1, []);
            end
            
            F_Total = F_Total - F_AddedMass;
            
        end
        
        function timeDomainSimReset (self)
            % reset the hydrodynamic system for transient simulation
            %
            % Syntax
            %
            %  timeDomainSimReset (hsys)
            %
            % Input
            %
            %  hsys - hydroSystem object
            %
            
            % call each body's timeDomainSimReset method
            for ind = 1:numel(self.hydroBodies)
                self.hydroBodies(ind).timeDomainSimReset ();
            end
            
            self.odeSimInitialised = false;
            
        end

        function n = get.nHydroBodies (self)
            
            n = numel (self.hydroBodies);

        end
        
        
    end
    
    
end