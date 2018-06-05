classdef linearPMMachinePTO < wsim.powerTakeOff
    % power take-off from relative linear displacement of two nodes
   
    properties (GetAccess = public, SetAccess = private)
        
        initialInternalTimeStep;
        maxInternalTimeStep;
        design;
        simoptions;
        
    end
    
    properties (GetAccess = private, SetAccess = private)
        
        mbdynForceObj;
        machineODESolver;
        initialPhaseCurrents; % Initial machine phase currents
        
    end
    
    methods
        
        function self = linearPMMachinePTO (reference_node, other_node, axisNum, design, simoptions, varargin)
            % construct a wsim.linearPMMachinePTO object
            %
            %
            % Syntax
            %
            % lpto = wsim.linearPMMachinePTO (reference_node, other_node, axisNum)
            % lpto = wsim.linearPMMachinePTO (..., 'Parameter', value)
            %
            % Description
            %
            % wsim.linearPMMachinePTO is a class representing a linear
            % power-take-off mechanism in a wave energy converter. It
            % facilitates sending the correct forces to an MBDyn multibody
            % simulation. wsim.linearPMMachinePTO applies forces between
            % two MBDyn nodes based on their relative displacement. Forces
            % are applied based on the relative displacement and velocity
            % of the two nodes along axis 3 in the reference frame of the
            % first node. It is assumed that the nodes motion is
            % constrained appropriately by other MBDyn elements (e.g. 
            %
            % Input
            %
            %  reference_node - mbdyn.pre.structuralNode6dof object
            %
            %  other_node - mbdyn.pre.structuralNode6dof object
            %
            %  axisNum - axis in the frame of the reference node. Forces
            %   will be applied to the node in a direction parallel to this
            %   axis.
            %
            %  design - 
            %
            %  simoptions - 
            %
            % Additional options my be supplied as parameter-value pairs.
            % The avaialable options are:
            %
            %  'InitialDisplacementZero' - optional true/false flag
            %    indicating whether the intial relative displacement (along
            %    axis 3 of the reference node) in the global frame should
            %    be taken as the reference point for displacement during
            %    the simulation, i.e. the PTO starts with an initial
            %    displacement of zero for the purposes of force calulation,
            %    and future displacement is measured relative to this
            %    initial position. If false, the raw position is used
            %    instead. Default is true if not supplied.
            %
            % Output
            %
            %  lpto - a wsim.linearPMMachinePTO
            %
            %
            %
            % See Also: wsim.powerTakeOff, wsim.linearPowerTakeOff
            %

            options.InitialDisplacementZero = true;
            options.ForceFcn = [];
            options.InitialPhaseCurrents = [0,0,0];
            options.InitialStep = [];
            options.MaxStep = [];
            
            options = parse_pv_pairs (options, varargin);
            
            if ~isempty (options.InitialStep)
                check.isNumericScalar (options.InitialStep, true, 'InitialStep', 1);
                assert (options.InitialStep > 0, 'InitialStep must be greater than zero');
            end
            if ~isempty (options.MaxStep)
                check.isNumericScalar (options.InitialStep, true, 'MaxStep', 1);
                assert (options.InitialStep > 0, 'MaxStep must be greater than zero');
            end
            
            lginfo.AvailableNames = { 'PTOTime', ...
                                      'Fpto', ...
                                      'EMF', ...
                                      'PhaseCurrents', ...
                                      'FaddE', ...
                                      'RelativeDisplacement', ...
                                      'RelativeVelocity' ...
                                    };
                              
            lginfo.IndepVars = { '', ...
                                 1, ...
                                 1, ...
                                 1, ...
                                 1, ...
                                 'Time', ...
                                 'Time', ...
                                 };
                              
            lginfo.Sizes = { [1,1], [1,1], [1,3], [1,3], [1,1], [1,1], [1,1] };
            
            
            lginfo.AxisLabels = { 'Time [s]', ...
                                  'Force [N]', ...
                                  'EMF [V]', ...
                                  'Phase Currents [A]', ...
                                  'Force [N]', ...
                                  'Relative Displacement [m]', ...
                                  'Relative Velocity [ms^{-1}]', ...
                                };
                            
            lginfo.Descriptions = repmat ({''}, size (lginfo.AxisLabels));
            
            lginfo.NAvailable = numel(lginfo.AvailableNames);
            
            self = self@wsim.powerTakeOff (reference_node, other_node, lginfo);
            
            self.mbdynForceObj = mbdyn.mint.twoNodeTranslationalForce ( ...
                                    reference_node, other_node, axisNum, ...
                                    'InitialDisplacementZero', options.InitialDisplacementZero, ...
                                    'ForceFcn', [] );
            
            self.initialPhaseCurrents = options.InitialPhaseCurrents;
            self.design = design;
            self.simoptions = simoptions;
            self.initialInternalTimeStep = options.InitialStep;
            self.maxInternalTimeStep = options.MaxStep;
            
        end
        
        function [FM, ptoforce, reldisp, relvel] = forceAndMoment (self, time)
            
            if time == self.simulationInfo.TStart
                
                FM = zeros (6,2);
                
            else
                
                [reldisp, relvel] = displacements (self.mbdynForceObj);

                solve (self.machineODESolver, time, [reldisp; relvel]);

                results = self.machineODESolver.getOutputs ();

                F = self.mbdynForceObj.force (results.Fpto);

                % need to add zero moments to forces
                FM = [ F; 
                       zeros(size (F)) ];
                
            end
            
        end
        
        function start (self, siminfo)
            % initialise the pto simulation
            
            start@wsim.powerTakeOff (self, siminfo);
            
            [reldisp, relvel] = self.mbdynForceObj.displacements ();
            
            innerodeevfcn = @(t, y, mc, flag) nestedodeforcefcn_linear (t, y, mc, flag, self.design, self.simoptions);

            % initial displacement and velocity of the generator
            interpdat = [reldisp; relvel];
            
            if isempty (self.initialInternalTimeStep)
                self.initialInternalTimeStep = self.simulationInfo.TStep/3; 
            end
                
            odeoptions = odeset ('InitialStep', self.initialInternalTimeStep);
            
            if ~isempty (self.maxInternalTimeStep)
                odeoptions = odeset (odeoptions, 'MaxStep', self.maxInternalTimeStep);
            end

            self.machineODESolver = ode.odesolver ( ...
                            innerodeevfcn, ...
                            self.simulationInfo.TStart, ...
                            self.initialPhaseCurrents, ...
                            interpdat, ...
                            odeoptions, ...
                            'Solver', 'ode15s', ...
                            'SaveSolutions', false ...
                                ... 'SplitFcn', @(flag, results, sol, mc, evalfcn) nestedsysresults_linear (flag, results, sol, mc, evalfcn, design, simoptions) ...
                                          );
            
        end
        
        function advanceStep (self, time, varargin)
            % advance to the next simulation time step
            %
            % Syntax
            %
            % advanceStep (pto)
            %
            % Description
            %
            % wsim.linearPMMachinePTO.advanceStep is intended to be called
            % when the simulation is ready to advance to the next time
            % step. It is always called when the linearPMMachinePTO is used
            % with the wsim.wecSim class to manage a simulation.
            % wsim.linearPMMachinePTO.advanceStep log data to the
            % wsim.logger object (internally it calls
            % wsim.powerTakeOff.logData). It also recalculates the last
            % values of the mahine simualtion at the current time step and
            % advances the ODE simulation of the PM machine to the next
            % time step, accepting the last values in to the simulation
            % history.
            %
            % Input
            %
            %  pto - wsim.powerTakeOff object
            %
            %
            % See Also: wsim.powerTakeOff.logData
            %
            
            [reldisp, relvel] = displacements (self.mbdynForceObj);
            
            if time == self.simulationInfo.TStart
                
%                 self.machineODESolver.acceptState ([reldisp; relvel]);
                
                % get the full output from the last generator sim
                self.internalVariables = self.machineODESolver.getOutputs ({}, false);
                
                self.internalVariables.PTOTime = time;
                
                self.internalVariables.PhaseCurrents = self.machineODESolver.sol.y';
                
            else
                % get the full output from the last generator sim
                self.internalVariables = self.machineODESolver.getOutputs ({}, true);
                
                % strip the first entry as after the first step it is a
                % duplicate of the final step of the previous simulation block
                fnames = fieldnames (self.internalVariables);
                for ind = 1:numel (fnames)
                    self.internalVariables.(fnames{ind})(1,:) = [];
                end
                
                self.internalVariables.PTOTime = self.machineODESolver.sol.x(2:end)';
            
                self.internalVariables.PhaseCurrents = self.machineODESolver.sol.y(:,2:end)';
                
                self.machineODESolver.acceptState ([reldisp; relvel]);
            
            end
            
            self.internalVariables.RelativeDisplacement = reldisp;
            
            self.internalVariables.RelativeVelocity = relvel;
            
            
            
            self.logData ();
            
        end
        
        function logData (self)
            % appends the internal variable data to the log
            %
            % Syntax
            %
            % logData (pto)
            %
            % Description
            %
            % logData appends the last recorded values of the internal
            % variables to a logger object.
            %
            % Input
            %
            %  pto - wsim.linearPMMachinePTO object
            %
            % See Also: wsim.powerTakeOff.loggingSetup,
            %           wsim.linearPMMachinePTO.advanceStep
            %
            
            if self.loggerReady
                for ind = 1:numel(self.loggingInfo.LoggedVarInds)
                    
                    fieldname = self.loggingInfo.AvailableNames{self.loggingInfo.LoggedVarInds(ind)};
                    
                    for indii = 1:size ( self.internalVariables.(fieldname), 1)
                        self.logger.logVal ( self.uniqueLoggingNames{self.loggingInfo.LoggedVarInds(ind)}, ...
                                             self.internalVariables.(fieldname)(indii,:) ...
                                           );
                    end
                end
            else
                error ('You have called logData, but logging has not been set up, have you called loggingSetup yet?');
            end
            
        end
        
    end
    
end