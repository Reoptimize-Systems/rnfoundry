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
%    large currents due to inductances. By default this ramp lasts for 2% of 
%    the specified PoleCount, but this can be modified by setting the 
%    RampPoles field in the simoptions structure.
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

        simoptions = setfieldifabsent (simoptions, 'RampPoles', ceil (simoptions.PoleCount * 2/100));

        simoptions = simsetup_ROTARY(design, simoptions.simfun, simoptions.finfun, ...
                                'RPM', simoptions.RPM, ...
                                'PoleCount', simoptions.PoleCount, ...
                                'RampPoles', simoptions.RampPoles, ...
                                'odeevfun', 'prescribedmotode_linear', ...
                                'simoptions', simoptions);

    end

%     simoptions.abstol = repmat(0.001, 1, design.Phases);

end