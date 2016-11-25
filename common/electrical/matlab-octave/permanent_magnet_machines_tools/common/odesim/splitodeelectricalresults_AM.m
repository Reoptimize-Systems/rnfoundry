function results = splitodeelectricalresults_AM (design, simoptions, flag, results, sol, odeinternals, phasecurrentcol, powermult, longresults)
% determines and stores useful intermediate electrical results during a
% split ode simulation of an electrical machine

    if nargin < 9
        longresults = false;
    end
    
    if flag == 0
        
        results.coilIsquaredsum = zeros(1,design.Phases);
        
        results.phaseIsquaredsum = zeros(1,design.Phases);

        results.coilEMFsquaredsum = zeros(1,design.Phases);
        
        results.phaseEMFsquaredsum = zeros(1,design.Phases);

        results.ICoilPeak = zeros(1,design.Phases);

        results.EMFPhasePeak = zeros(1,design.Phases);

        results.EnergyLoadTotal = zeros(1,design.Phases);
        
        results.EnergyPhaseRTotal = zeros(1,1);

        results.GridPowersum = zeros(1,design.Phases);
        % the peak total exported power
        results.PowerLoadPeak = 0;
        
        results.maxT = [];
        
        return;

    elseif flag == 1
    
        if nargin < 5
            % this factor is used to multiply the power by, say, the number of
            % stages in the machine
            powermult = 1;
        end

        % Determine some interesting design outputs common to all machines, all
        % results here are on a per-generator basis

        % extract the first coil current from the solution, the current in
        % sol.y(simoptions.ODEPhaseCurrentCol,:) is the phase current of the
        % first phase, so we divide by the number of parallel branches to
        % obtain the current in the coils
        phaseCurrent = sol.y(phasecurrentcol,:)';

        % calculate the coil current
        coilCurrent = phaseCurrent ./ design.Branches;

        % Get the square of the current values
        coilCurrentSquared = realpow(coilCurrent, 2);
        phaseCurrentSquared = realpow(phaseCurrent, 2);

        % The integral of the square of the current, can be used to
        % determine the rms current later
        results.coilIsquaredsum = results.coilIsquaredsum + trapz(sol.x(:), coilCurrentSquared);
        results.phaseIsquaredsum = results.phaseIsquaredsum + trapz(sol.x(:), phaseCurrentSquared);

        % The integral of the square of the EMF, can be used to determine the
        % rms EMF later
        results.coilEMFsquaredsum = results.coilEMFsquaredsum + trapz(sol.x(:), realpow(odeinternals.EMF, 2));

        phaseEMF = odeinternals.EMF .* design.CoilsPerBranch;
        results.phaseEMFsquaredsum = results.phaseEMFsquaredsum + trapz(sol.x(:), realpow(phaseEMF, 2));

        % peak branch (coil current) (phase current divided by number of
        % parallel branches))
        results.ICoilPeak = max([results.ICoilPeak; abs(coilCurrent)], [], 1);

        % peak phase EMF (coil EMF times the number of coils per branch)
        results.EMFPhasePeak = max([results.EMFPhasePeak; abs(phaseEMF)], [], 1);

        % Calculate the power output at each time step
        GridPower = phaseCurrentSquared .* design.LoadResistance * powermult;

        % The integral of the energy dissipated in the load resistance,
        % i.e. the integral of the grid power
        results.EnergyLoadTotal = results.EnergyLoadTotal + trapz(sol.x(:), GridPower);

        % The integral of the energy dissipated in the phase resistances
        PhasePower = sum(realpow(phaseCurrent,2) .* odeinternals.RPhase, 2);
        results.EnergyPhaseRTotal = results.EnergyPhaseRTotal + trapz(sol.x(:), PhasePower);

        % Determine the peak power output
        results.PowerLoadPeak = max(abs([results.PowerLoadPeak; sum(GridPower,2)]));

        % get the maximum time of the sim
        results.maxT = sol.x(end);
        
        if longresults
            if results.block == 1
%                 results.Time = sol.x';
                results.EMF = phaseEMF;
                results.PhaseCurrents = sol.y(simoptions.ODESim.NestedSim.SolutionComponents.PhaseCurrents.SolutionIndices,:).';
            else
%                 results.Time = [results.Time; sol.x(2:end)'];
                results.EMF = [results.EMF; phaseEMF(2:end,:)];
                results.PhaseCurrents = [results.PhaseCurrents;
                    sol.y(simoptions.ODESim.NestedSim.SolutionComponents.PhaseCurrents.SolutionIndices,2:end).'];
            end
        end
        
    elseif flag == 2
        
%         % determine the actual simulation time
%         design.SimTimeSpan = results.maxT - simoptions.ODESim.TimeSpan(1);
% 
%         % Get the peak coil current and current densities
%         [design.ICoilPeak, maxIind] = max(results.ICoilPeak);
% 
%         design.JCoilPeak = design.ICoilPeak / design.ConductorArea;
% 
%         % Calculate the rms coil current from the supplied current integral
%         design.ICoilRms = sqrt(results.coilIsquaredsum ./ ( design.SimTimeSpan ));
% 
%         % The rms coil current density
%         design.JCoilRms = design.ICoilRms(maxIind) / design.ConductorArea;
% 
%         % The max current density in the conductor
%     %     results.JCoilPeak = results.ICoilPeak(maxIind) /
%     %     design.ConductorArea;
% 
%         % Calculate the rms phase EMF from the supplied EMF integral
%         design.EMFPhaseRms = sqrt(results.phaseEMFsquaredsum(maxIind) ./ ( design.SimTimeSpan ));
% 
%         % Get the peak phase current 
%         design.IPhasePeak = design.ICoilPeak .* design.Branches;
% 
%         % Get the rms phase current
%         design.IPhaseRms = sqrt(results.phaseIsquaredsum(maxIind)  ./ ( design.SimTimeSpan ));
% 
%         % Store the peak phase EMF in the design structure
%         design.EMFPhasePeak = max(results.EMFPhasePeak);
% 
%         % Store the total extracted energy summed for all generators in the
%         % design structure
%         design.EnergyLoadTotal = sum(results.EnergyLoadTotal) * simoptions.NoOfMachines;
% 
%         design.EnergyPhaseRTotal = results.EnergyPhaseRTotal * simoptions.NoOfMachines;
% 
%         % Calculate the mean power output from the supplied power integral
%         % of the simulation
%         design.PowerLoadMean = design.EnergyLoadTotal / (design.SimTimeSpan);
% 
%         % put the peak recorded power in the design structure
%         design.PowerLoadPeak = results.PowerLoadPeak * simoptions.NoOfMachines;
% 
%         % store the mean power losses in the phase resistances
%         design.PowerPhaseRMean = design.EnergyPhaseRTotal / (design.SimTimeSpan);
% 
%         design.PowerSystemMean = design.PowerPhaseRMean + design.PowerLoadMean;
    
    end

end