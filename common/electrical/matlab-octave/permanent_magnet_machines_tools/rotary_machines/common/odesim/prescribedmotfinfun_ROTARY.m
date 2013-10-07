function [design, simoptions] = prescribedmotfinfun_ROTARY(design, simoptions, finfun)
% performs postprocessing of a machine design in preparation for a rotary
% proescribed motion simulation
%
% Syntax
%
% [design, simoptions] = prescribedmotfinfun_ROTARY(design, simoptions, finfun)
%
% Input
%
%   design and simoptions are structures containing the particulars of the
%    design under consideration. For prescribedmotfinfun_ROTARY the
%    simoptions structure may contain the fields 'PoleCount' and 'RPM'. If
%    present a set of times, positions and velocities will be added to the
%    simoptions structure which yield a constant velocity simulation
%    covering the displacement equivalent to the number of Poles in
%    polecount at the RPM in the 'RPM' field. A linear ramp up in speed from
%    zero to the specified RPM is prepended to this simulation to avoid
%    large currents due to inductances.
%
%   finfun is a function handle or string containing a function name which
%     is the preprocessing function for the design. 
%
% Output
%
% 

    if ~all(isfield(design, {'slm_psidot'}))
        % In this case we assume we have not already run the finalisation
        % code on this design and must do so
        [design, simoptions] = feval(finfun, design, simoptions);
    end
    
    if all(isfield(simoptions, {'PoleCount', 'RPM'}))
        
        simoptions = simsetup_ROTARY(design, simoptions.simfun, simoptions.finfun, ...
                                'RPM', simoptions.RPM, ...
                                'PoleCount', simoptions.PoleCount, ...
                                'odeevfun', 'prescribedmotode_linear', ...
                                'simoptions', simoptions);
                          
        
        % add a linear speed ramp up over 5 Poles to reduce the starting
        % currents due to inductance
        nramppoles = 5;
        rampa = simoptions.omegaT(1)^2 / (2 * nramppoles * design.PoleWidth);
        rampTmax = simoptions.omegaT(1) / rampa;
        rampT = linspace(0, rampTmax, 15);
        rampomegaT = rampa .* rampT;
        rampthetaT = 0.5 * rampa .* rampT.^2;
        
        simoptions.omegaT = [rampomegaT(1:end-1), simoptions.omegaT];
        simoptions.thetaT = [rampthetaT(1:end-1), simoptions.thetaT + rampthetaT(end)];
        simoptions.drivetimes = [rampT(1:end-1), simoptions.drivetimes + rampT(end)];
        
        simoptions.tspan = simoptions.drivetimes([1, end]);

        simoptions.maxstep = (simoptions.tspan(2) - simoptions.tspan(1)) / (40 * simoptions.PoleCount);
        
        % construct a piecewise polynomial interpolation of the position
        % and velocity data
        simoptions.pp_thetaT = interp1(simoptions.drivetimes,simoptions.thetaT,'cubic','pp');
        simoptions.pp_omegaT = interp1(simoptions.drivetimes,simoptions.omegaT,'cubic','pp');
        
    end
    
%     simoptions.abstol = repmat(0.001, 1, design.Phases);
        
end