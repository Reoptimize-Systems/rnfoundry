function [results] = splitodeelectricalresults_AM(flag, design, simoptions, results, sol, odeinternals, powermult)
% determines and stores useful intermediate electrical results during a
% split ode simulation of an electrical machine

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
        
    else
        if nargin < 5
            powermult = 1;
        end
    end

    % Determine some interesting design outputs common to all machines, all
    % results here are on a per-generator basis

    % extract the first coil current from the solution, the current in
    % sol.y(simoptions.ODEPhaseCurrentCol,:) is the phase current of the
    % first phase, so we divide by the number of parallel branched to
    % obtain the current in the coils
    phaseCurrent = sol.y(simoptions.ODEPhaseCurrentCol:simoptions.ODEPhaseCurrentCol-1+design.Phases,:)';
    
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

end