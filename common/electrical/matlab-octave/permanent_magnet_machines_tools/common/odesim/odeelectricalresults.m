function design = odeelectricalresults(T, Iphase, EMF, RPhase, design, simoptions)
% odeelectricalresults: calculates the electrical outputs, and exported
% power from an ode simulation of a generator

    % Determine some interesting electrical design outputs
    
    % we should use the phase that produced the highest current in some
    % calculations
    [design.IPhasePeak,maxIind] = max(max(abs(Iphase), [], 1));
    
    % rms phase current, use phase current which has the max value
    design.IPhaseRms = contrms(T, Iphase(:,maxIind));
    
    % rms branch (coil) current
    design.ICoilRms = contrms(T, Iphase(:,maxIind) ./ design.Branches);
    design.JCoilRms = design.ICoilRms  / design.ConductorArea;
    
    % rms phase EMF (coil EMF times the number of coils per branch)
    design.EMFPhaseRms = contrms(T, EMF(:,maxIind) .* design.CoilsPerBranch);
    
    % peak branch (coil current) (phase current divided by number of
    % parallel branches))
    design.ICoilPeak = design.IPhasePeak ./ design.Branches;
    design.JCoilPeak = design.ICoilPeak / design.ConductorArea;
    
    % peak phase EMF (coil EMF time the number of coils per branch)
    design.EMFPhasePeak = max(abs(EMF(:))) * design.CoilsPerBranch;

    if numel(T) > 1
        
        % calculate the power from the phase currents in the load
        loadPower = sum(realpow(Iphase,2), 2) * design.LoadResistance ...
                        * design.NStages * simoptions.NoOfMachines;
        
        design.EnergyLoadTotal = sum(trapz(T, loadPower));
        
        design.PowerLoadMean = contmean(T, loadPower);
        
        phasePower = sum(bsxfun (@times, realpow(Iphase,2), RPhase), 2) ...
                        .* design.NStages .* simoptions.NoOfMachines;
        
        design.EnergyPhaseRTotal = trapz(T, phasePower);

        design.PowerPhaseRMean = contmean(T, phasePower);

        design.PowerSystemMean = design.PowerPhaseRMean + design.PowerLoadMean;
        
        design.PowerLoadPeak = max(loadPower);
        
        design.PowerFactorEstimate = (design.Phases ...
                                        .* design.IPhaseRms.^2 ...
                                        .* (design.PhaseResistance(end)+design.LoadResistance)) ...
                                      ./ (design.Phases * design.EMFPhaseRms * design.IPhaseRms);

%         
%         Y = fft(y,251);
%         Pyy = Y.*conj(Y)/251;
%         f = 1000/251*(0:127

%         % use Lomb normalized periodogram to find the main frequency,
%         % cannot use fft (without interpolating) as data is not uniformly
%         % spaced
%         [design.Periodogram,design.PeriodogramFrequencies] = fastlomb(EMF,T,0,2);
%         
%         % find the main harmonic frequency and store it
%         [junk,I] = max(design.Periodogram);
%         
%         design.ElectricalFrequency = design.PeriodogramFrequencies(I);
        
        % use this frequency to calculate the power factor
        % design.PowerFactor = rlcpowerfactor(design.R(1,1), design.L(1,1), 0, 2 * pi * design.ElectricalFrequency);
        
    else
        
        design.EnergyLoadTotal = 1e-4;
        design.EnergyLoadMean = 1e-4;
        design.PowerLoadMean = 1e-4;
        design.PowerLoadPeak = 1e-4;
        design.PowerSystemMean  = 1e-4;
        design.PowerPhaseRMean = 1e-4;
        
    end

end