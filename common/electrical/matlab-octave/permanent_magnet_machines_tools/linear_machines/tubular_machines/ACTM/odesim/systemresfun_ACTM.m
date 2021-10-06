function [results, design] = systemresfun_ACTM(T, Y, design, simoptions)

    % Now obtain internally calculated values by recalling the function
    % with the results, first preallocating the arrays. This is necessary
    % as the ode solver used may take steps while choosing step sizes which
    % do not form part of the solution      
    
    % We may not want to recalculate and plot every single solution step,
    % so we use the skip value to allow us to skip every x solution points,
    % e.g. skip = 2 would calculate only every other solution point.
    simoptions.ODESim.ResultsTSkip = 1;
    
    % Now preallocate arrays of the correct sizes
    results.ydot = zeros(ceil(size(Y,1)/simoptions.ODESim.ResultsTSkip), size(simoptions.ODESim.InitialConditions,2));
    results.dpsidxR = zeros(ceil(length(T)/simoptions.ODESim.ResultsTSkip), 3);
    results.EMF = zeros(ceil(length(T)/simoptions.ODESim.ResultsTSkip), 3);
    results.Ffea = zeros(ceil(length(T)/simoptions.ODESim.ResultsTSkip), 1);
    results.xT = zeros(ceil(length(T)/simoptions.ODESim.ResultsTSkip), 1);
    results.vT = zeros(ceil(length(T)/simoptions.ODESim.ResultsTSkip), 1);
    results.excitation_force_heave = zeros(ceil(length(T)/simoptions.ODESim.ResultsTSkip), 1);
    results.excitation_force_surge = zeros(ceil(length(T)/simoptions.ODESim.ResultsTSkip), 1);
    results.radiation_force_heave = zeros(ceil(length(T)/simoptions.ODESim.ResultsTSkip), 1);
    results.radiation_force_surge = zeros(ceil(length(T)/simoptions.ODESim.ResultsTSkip), 1);
    results.buoyancy_force = zeros(ceil(length(T)/simoptions.ODESim.ResultsTSkip), 1);
    results.FBDh = zeros(ceil(length(T)/simoptions.ODESim.ResultsTSkip), 1);
    results.FBDs = zeros(ceil(length(T)/simoptions.ODESim.ResultsTSkip), 1);
    results.Ffea_heave = zeros(ceil(length(T)/simoptions.ODESim.ResultsTSkip), 1);
    results.Ffea_surge = zeros(ceil(length(T)/simoptions.ODESim.ResultsTSkip), 1);
    % Recalculate the final solution values
    
    k = 0;
    for z = 1:simoptions.ODESim.ResultsTSkip:length(T)

        k = k + 1;

        [results.ydot(k,:),...
            results.EMF(k,:),...
            results.Ffea(k,1),...
            results.xT(k,1),...
            results.vT(k,1),...
            results.Ffea_heave(k,1),...
            results.Ffea_surge(k,1),...
            results.FBDh(k,1),...
            results.FBDs(k,1),...
            results.bouyancy_force(k,1),...
            results.excitation_force_heave(k,1),...
            results.excitation_force_surge(k,1),...
            results.radiation_force_heave(k,1),...
            results.radiation_force_surge(k,1),...
            results.dpsidxR(k,:)] = systemode_ACTM(T(z,1), Y(z,:)', design, simoptions);
        
    end

    % Calculate the wave heights
    results.wave_height = zeros(1, length(T));

    for k = 1:length(T)
        time_vector = T(k) * ones(1, length(simoptions.BuoySim.SeaParameters.phase));
        results.wave_height(k) = sum(real(simoptions.BuoySim.SeaParameters.amp .* exp(-i .* (simoptions.BuoySim.SeaParameters.sigma .* time_vector - simoptions.BuoySim.SeaParameters.phase))));
    end
    
    design.xAmax = 0;

    % Actual armature length that would be required to achieve the
    % output, this is not used in the cost calculations
    peakxT = max(abs(results.xT));

    design.minArmLength = 2 * peakxT + (design.Poles(1) * design.Wp);

    design.minArmPoles = ceil(design.minArmLength ./ design.Wp);

    design.minArmLength = design.minArmPoles * design.Wp;

    % Determine some interesting design outputs
    design.Irms = rms(interp1(T, Y(:,5), 0:max(T)/(length(T)*2):max(T)));
    design.Jrms = design.Irms  / design.conductorArea;
    design.EMFrms = rms(interp1(T, results.EMF(:,1), 0:max(T)/(length(T)*2):max(T)));
    design.Ipeak = max(abs(Y(:,5)));
    design.EMFpeak = max(abs(results.EMF(:,1)));
    design.maxJ = max(abs(Y(:,5))) / design.conductorArea;
    design.TotalGridEnergy = trapz(T, Y(:,5).^2 .* design.LoadResistance) * design.Poles(1) * 3;
    design.AverageEnergy = design.TotalGridEnergy ./ (max(T) - min(T));
    design.GridMeanPower = mean(interp1(T, Y(:,5), 0:max(T)/(length(T)*2):max(T)).^2 * design.LoadResistance) * design.Poles(1) * 3;
    design.peakPower = max(Y(:,5).^2 .* design.LoadResistance) * design.Poles(1) * 3;

end
