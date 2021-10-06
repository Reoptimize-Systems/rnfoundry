function [design, simoptions] = systemfinfun_linear(design, simoptions, finfun)
% systemfinfun_linear: finalises the machine design and buoy simulation
% setup for a linear machine and heaving buoy simulation

    if ~all(isfield(design, {'slm_fluxlinkage'}))
        % In this case we assume we have not already run the finalisation
        % code on this design and must do so
        [design, simoptions] = feval (finfun, design, simoptions);
    end
    
    % set the maximum allowed displacement of the translator, the default
    % is infinite. If this is exceeded end stops are used to stop the
    % motion
    simoptions.BuoySim = setfieldifabsent (simoptions.BuoySim, 'maxAllowedxT', inf);
    
    simoptions.BuoySim = setfieldifabsent (simoptions.BuoySim, 'buoy', []);
    
    % Set up buoy simulation
    simoptions.BuoySim = buoysimsetup (simoptions.BuoySim.buoy, simoptions.BuoySim);
    
    % copy over the buoy ode simulation components info created by
    % buoysimsetup
    fnames = fieldnames (simoptions.BuoySim.ODESim.SolutionComponents);
    for ind = 1:numel (fnames)
        simoptions.ODESim.SolutionComponents.(fnames{ind}) ...
            = simoptions.BuoySim.ODESim.SolutionComponents.(fnames{ind});
    end
    
    % create an empty Buoy structure if not present
    design = setfieldifabsent (design, 'Buoy', struct ());
    
    % copy over the buoy directory used in the simulation
    design.Buoy.Designation = simoptions.BuoySim.buoy;
    
    % add zero WEC friction if not supplied
    design.Buoy = setfieldifabsent (design.Buoy, 'WECFriction', 0);

end