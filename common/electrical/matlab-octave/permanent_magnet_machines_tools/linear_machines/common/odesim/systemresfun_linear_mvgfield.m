function [results, design] = systemresfun_linear_mvgfield(T, Y, design, simoptions)

    % Now obtain internally calculated values by recalling the function
    % with the results, first preallocating the arrays. This is necessary
    % as the ode solver used may take steps while choosing step sizes which
    % do not form part of the solution

    % We may not want to recalculate and plot every single solution step,
    % so we use the skip value to allow us to skip every x solution points,
    % e.g. skip = 2 would calculate only every other solution point.
    simoptions.skip = 1;

    if ~isfield(simoptions, 'ODEPhaseCurrentCol')
        simoptions.ODEPhaseCurrentCol = 7;
    end
    
    [results, design] = resfun_linear(T, Y, design, simoptions);
    
%     % Calculate the wave heights
%     results.wave_height = zeros(1, length(T));
%
%     for k = 1:length(T)
%         time_vector = T(k) * ones(1, length(simoptions.BuoySim.SeaParameters.phase));
%         results.wave_height(k) = sum(real(simoptions.BuoySim.SeaParameters.amp .* exp(-i .* (simoptions.BuoySim.SeaParameters.sigma .* time_vector - simoptions.BuoySim.SeaParameters.phase))));
%     end
    
    % FF = [Fs, FLinearDrag, FdragF, FEAFy];
    
    results.Fs = results.FaddF(:,1);
    
    results.Fsnap = results.FaddF(:,end);
    
%     if numel(results.FaddF) > 2
%         results.FfF = results.FaddF(:,2);
%     end
    
    if numel(results.FaddF) > 2
        results.FLinearDrag = results.FaddF(:,2);
    end
    
    if numel(results.FaddF) > 3
        results.FdragF = results.FaddF(:,3);
    end

    % Actual translator length that would be required to achieve the
    % output, this is not used in the cost calculations
    peakxT = max(results.xT);
    troughxT = min(results.xT);
    peakxF = max(Y(:,5));
    troughxF = min(Y(:,5));
    
%     design.xAmax = max(abs(peakxF, troughxF));

    design.minTransLength = 2 * max(peakxT - troughxF, peakxF - troughxT) + (design.Poles(1) * design.PoleWidth);
    %design.minTransLength = max(max(abs(Y(:,6) - results.xT)),max(results.xT)) + (design.Poles(1) * design.PoleWidth);
    design.minTransPoles = ceil(design.minTransLength ./ design.PoleWidth);

    design.minTransLength = design.minTransPoles * design.PoleWidth;

    design.extraTMass = 1.1 * max(abs(results.Fpto)) / simoptions.BuoySim.BuoyParameters.g;

%     design.minTransMass = transMass_Snapper(design);

%     design = odeelectricalresults(T, Y(:,7), results.EMF(:,1), design, simoptions);

    design.xFmax = max(abs(Y(:,5)));
    design.xFrms = rms(interp1(T, Y(:,5), 0:max(T)/(length(T)*2):max(T)));

%         simoptions.BuoySim.BuoyParameters.mass = simoptions.BuoySim.BuoyParameters.draft * pi * simoptions.BuoySim.BuoyParameters.a^2 * simoptions.BuoySim.BuoyParameters.rho;

end
        