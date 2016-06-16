function [results, design] = systemresfun_linear_mvgarm(T, Y, design, simoptions)

    % Now obtain internally calculated values by recalling the function
    % with the results, first preallocating the arrays. This is necessary
    % as the ode solver used may take steps while choosing step sizes which
    % do not form part of the solution

    % We may not want to recalculate and plot every single solution step,
    % so we use the skip value to allow us to skip every x solution points,
    % e.g. skip = 2 would calculate only every other solution point.
%     simoptions.skip = 1;

    if ~isfield(simoptions, 'ODEPhaseCurrentCol')
        simoptions.ODEPhaseCurrentCol = 7;
    end

    [results, design] = resfun_linear(T, Y, design, simoptions);

    % from forcefcn_linear_mvgarm: 
    % FaddEBD = [Fs, FdragA, FEAFy, FLinearDrag, FfT, FdragT, Fes, Floss]
    % spring force
    results.Fs = results.FaddEBD(:,1);
    % fluid drag forces on armature
    results.FdragA = results.FaddEBD(:,2);
    % magnetic coupling force
    results.Fsnap = results.FaddEBD(:,3);
    % artificial linear drag force 
    results.FLinearDrag = results.FaddEBD(:,4);
    % translator friction
    results.FfT = results.FaddEBD(:,5);
    % end stop forces
    results.Fes = results.FaddEBD(:,7);
    % forces due to losses in the machine
    results.FLoss = results.FaddEBD(:,8);
    
%     % Calculate the wave heights
%     results.wave_height = zeros(1, length(T));
% 
%     for k = 1:length(T)
%         time_vector = T(k) * ones(1, length(simoptions.BuoySim.SeaParameters.phase));
%         results.wave_height(k) = sum(real(simoptions.BuoySim.SeaParameters.amp .* exp(-i .* (simoptions.BuoySim.SeaParameters.sigma .* time_vector - simoptions.BuoySim.SeaParameters.phase))));
%     end

    design.xAmax = 0;

    % Actual translator length that would be required to achieve the
    % output, this is not used in the cost calculations
    peakxT = max(results.xT);
    troughxT = min(results.xT);
    peakxA = max(Y(:,5));
    troughxA = min(Y(:,5));

    design.minLongMemberLength = 2 * max(peakxT - troughxA, peakxA - troughxT) + (max(design.Poles) * design.PoleWidth);
    
    design.minLongMemberPoles = ceil(design.minLongMemberLength ./ design.PoleWidth);

    design.minLongMemberLength = design.minLongMemberPoles * design.PoleWidth;

    design.extraFptoMass = 1.1 * max(abs(results.Fpto)) / simoptions.BuoySim.BuoyParameters.g;

    design.xAmax = max(abs(Y(:,5)));
    design.xArms = contrms(T, Y(:,5));
    design.xAmean = contmean(T, Y(:,5));
    design.vRmax = max(abs(results.vT - Y(:,6)));

end
        