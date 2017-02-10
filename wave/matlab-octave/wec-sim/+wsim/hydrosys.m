classdef hydrosys < handle
    % class representing a collection of bodies interacting with water
    % waves and each other hydrodynamically
    %
    %
    
    properties (GetAccess = private, SetAccess = private)
        
        simu;        % simulation setup
        waves;       % wave setup
        hydroBodies; % array of hydrobody objects
        hydroBodyInds = []; % array of indices of the hydrobodies
        
        odeSimInitialised = false;
        
    end
    
    
    methods
        
        function self = hydrosys (waves, simu, hydrobodies)
            % hydrosys class constructor
            %
            % Syntax
            %
            % hsys = hydrosys (waves, simu)
            % hsys = hydrosys (..., hydrobodyfiles)
            % 
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
            %  hsys - hydrosys object
            % 
            
            
            self.simu = simu;
            self.waves = waves;
            
            self.odeSimInitialised = false;
            
            % Waves and Simu: check inputs
            self.waves.checkinputs ();
            self.simu.checkinputs ();
            
            % initialise the number of WEC bodies in the sim to 0, this
            % will be incremented as they are added
            self.simu.numWecBodies = 0;
            
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
            %  hsys - hydrosys object
            %
            %  hydrobodies - array of one or more hydrobody objects to be
            %    added to the system
            %
            
            
            if self.odeSimInitialised
                
                error ('An ode simulation has already been initialised so you cannot add new bodies. To add new bodies first call the odeSimReset method.')
                
            else
                
                if ~isa (hydrobodies, 'wsim.hydrobody')
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
                    
                    % add the new hydrobodies to the collection
                    self.hydroBodies = [ self.hydroBodies, hydrobodies(ind) ];
                
                    % increment the number of bodies in simu
                    self.simu.numWecBodies = self.simu.numWecBodies + 1;
                    
                    nextbodyind = nextbodyind + 1;

                end
            
            end
            
        end
        
        function initialiseHydrobodies (self, varargin)
            
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
                
                self.hydroBodies(bodyind).checkinputs ();
                
                % Determine if hydro data needs to be reloaded from h5
                % file, or if hydroData was stored in memory from a
                % previous run.
                if options.MultiConditionRun == true ...
                        && self.simu.reloadH5Data == 0 ...
                        && options.MultiConditionRunIndex > 1 
                    
                    self.hydroBodies(bodyind).loadHydroData (hydroData(bodyind));
                    
                else
                    
                    self.hydroBodies(bodyind).readH5File;
                    
                end
                
                self.hydroBodies(bodyind).bodyTotal = self.simu.numWecBodies;
                
                if self.simu.b2b == 1
                    self.hydroBodies(bodyind).lenJ = zeros (6*self.hydroBodies(bodyind).bodyTotal, 1);
                else
                    self.hydroBodies(bodyind).lenJ = zeros (6, 1);
                end
                
            end 

        end
        
        function odeSimSetup (self)
            % set up the hydrodynamic system in preparation for performaing
            % a transient ode solution
            %
            % Syntax
            %
            %  odeSimSetup (hsys)
            %
            % Input
            %
            %  hsys - hydrosys object
            
            % simulation setup
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
            if (self.simu.nlHydro >0) || (self.simu.paraview == 1)
                for bodyind = 1:length(self.hydroBodies(1,:))
                    self.hydroBodies(bodyind).bodyGeo (body(bodyind).geometryFile)
                end
            end
            
            % call each body's odeSimSetup method
            for bodyind = 1:numel(self.hydroBodies)
                
                self.hydroBodies(bodyind).odeSimSetup (self.waves, self.simu, self.hydroBodyInds(bodyind));
                
            end
            
            for bodyind = 1:numel(self.hydroBodies)
                self.hydroBodies(bodyind).adjustMassMatrix(self.simu.adjMassWeightFun,self.simu.b2b);
            end
            
            self.odeSimInitialised = true;
            
        end
        
        function [forces, out] = hydroForces (self, t, x, vel, accel, elv)
            
            forces = nan * ones (6, self.simu.numWecBodies );
            
            for bodyind = 1:numel(self.hydroBodies)
                [forces(:,bodyind), out(bodyind)] = hydroForces (self.hydroBodies(bodyind), ...
                                                                 t, ...
                                                                 x(:,bodyind), ...
                                                                 vel, ...
                                                                 accel, ...
                                                                 elv );
            end
            
        end
        
        function [F_Total, F_AddedMass] = correctAddedMassForce (self, forceTotal, forceAddedMass, accel)
            % recalcualte the added mass and total forces on bodies
            
            F_Total = forceTotal + forceAddedMass;
            F_AddedMass = nan * ones (size (forceAddedMass));
            
            for bodyind = 1:numel(self.hydroBodies)
                self.hydroBodies(bodyind).restoreMassMatrix ();
                F_AddedMass(:,:,bodyind) = self.hydroBodies(bodyind).forceAddedMass(accel(:,:,bodyind), self.simu.b2b);
            end
            
            F_Total = F_Total - F_AddedMass;
            
        end
        
%         function odeSimReset (self)
%             % reset the hydrodynamic system for transient simulation
%             %
%             % Syntax
%             %
%             %  odeSimSetup (hsys)
%             %
%             % Input
%             %
%             %  hsys - hydrosys object
%             %
%             
%             % call each body's odeSimReset method
%             for ind = 1:numel(self.hydroBodies)
%                 self.hydroBodies(ind).odeSimReset ();
%             end
%             
%             self.odeSimInitialised = false;
%             
%         end
        
        
    end
    
    
end