function [design, simoptions] = prescribedmotfinfun_linear(design, simoptions, finfun)

    [design, simoptions] = prescribedmotfinfun_AM (design, simoptions, finfun);
    
%     simoptions.ODEPhaseCurrentCol = 1;
%     
%     % now set up the linear simulation parameters based on the desired
%     % velocity
%     if isfield(simoptions, 'PoleCrossings') && isfield(simoptions, 'Velocity')
%         
%         simoptions.ODESim.TimeSpan = [0, simoptions.PoleCrossings * design.PoleWidth / simoptions.Velocity];
%         simoptions.xT = [0, simoptions.PoleCrossings*design.PoleWidth];
%         simoptions.vT = [simoptions.Velocity, simoptions.Velocity];
%         simoptions.drivetimes = [0, simoptions.xT(end) / simoptions.Velocity];
%         
%     elseif ~isfield(simoptions, 'xT') && isfield(simoptions, 'vT')
%         
%         % integrate the velocity profile to get the displacement profile
%         simoptions.xT = cumtrapz(simoptions.drivetimes, simoptions.vT);
%         
%     elseif ~isfield(simoptions, 'xT') && isfield(simoptions, 'vT')
%         
%         % fit an slm model to the data
%         xTslm = slmengine(simoptions.drivetimes, simoptions.xT, 'knots', max(5, round(numel(simoptions.xT)/5)));
%         % find the derivative to get the velocity profile
%         simoptions.vT = slmeval(simoptions.drivetimes, xTslm, 1, false);
%         
%     elseif ~all(isfield(simoptions, {'drivetimes','xT','vT'}))
%         error('Insufficient options provided to set up prescribed motion sim.')
%     elseif ~isfield(simoptions, 'tspan')
%         simoptions.ODESim.TimeSpan = [simoptions.drivetimes(1), simoptions.drivetimes(end)];
%     end
%     
%     % construct a piecewise polynomial interpolation of the position
%     % and velocity data
%     simoptions.pp_xT = interp1(simoptions.drivetimes,simoptions.xT,'cubic','pp');
%     simoptions.pp_vT = interp1(simoptions.drivetimes,simoptions.vT,'cubic','pp');
%     
%     simoptions = setfieldifabsent(simoptions, 'IC', zeros(1, design.Phases));

end