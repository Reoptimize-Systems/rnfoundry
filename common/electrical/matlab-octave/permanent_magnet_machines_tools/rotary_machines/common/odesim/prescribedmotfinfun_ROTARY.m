function [design, simoptions] = prescribedmotfinfun_ROTARY(design, simoptions, finfun)
% performs postprocessing of a machine design in preparation for a rotary
% prescribed motion ode simulation
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

    [design, simoptions] = prescribedmotfinfun_AM (design, simoptions, finfun);
    
    % don't do rotary sim setup unless requested
    simoptions = setfieldifabsent (simoptions, 'DoRotarySimSetup', false);
    
    if simoptions.DoRotarySimSetup
        
        simoptions = setfieldifabsent (simoptions, 'RampPoles', ceil (simoptions.PoleCount * 2/100));
        simoptions.ODESim = setfieldifabsent (simoptions.ODESim, 'ForceAddPhaseCurrentODESolutionComps', true);

        simoptions = simsetup_ROTARY(design, simoptions.ODESim.PreProcFcn, simoptions.ODESim.PostPreProcFcn, ...
                                'RPM', simoptions.RPM, ...
                                'PoleCount', simoptions.PoleCount, ...
                                'RampPoles', simoptions.RampPoles, ...
                                'EvalFcn', simoptions.ODESim.PreProcFcn, ...
                                'ForceAddPhaseCurrentODESolutionComps', simoptions.ODESim.ForceAddPhaseCurrentODESolutionComps, ...
                                'simoptions', simoptions);

    end
    
    % set up power converter properties, if there is one and the required
    % data is also available
    if isfield (design, 'MachineSidePowerConverter') ...
         && isfield (simoptions, 'omegaT')
     
        maxomega = max (abs (simoptions.omegaT));
        
        % the minum DC link voltage is the natural rectification voltage at
        % the maximum expected speed
        minVdc = (3 * sqrt(6) / pi) * ...
                 peakemfest_ROTARY (design.FluxLinkagePhasePeak, maxomega, design.Poles/2) / sqrt(2);
             
        if ~isfield (design.MachineSidePowerConverter, 'Vdc')
            % make the DC link voltage 10% bigger than the minimum possible
            design.MachineSidePowerConverter.Vdc = 1.1 * minVdc;
            
        elseif design.MachineSidePowerConverter.Vdc < minVdc
            
            warning ('Machine side power converter DC voltage is lower than the minimum reccomended value of %gV', minVdc)
            
        end
        
    end

end