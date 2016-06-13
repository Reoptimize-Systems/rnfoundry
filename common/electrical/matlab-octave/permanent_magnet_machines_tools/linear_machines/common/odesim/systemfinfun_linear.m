function [design, simoptions] = systemfinfun_linear(design, simoptions, finfun)
% systemfinfun_linear: finalises the machine design and buoy simulation
% setup for a linear machine and heaving buoy simulation

    if ~all(isfield(design, {'slm_fluxlinkage'}))
        % In this case we assume we have not already run the finalisation
        % code on this design and must do so
        [design, simoptions] = feval(finfun, design, simoptions);
    end
    
    % set the maximum allowed displacement of the translator, the default
    % is infinite. If this is exceeded end stops are used to stop the
    % motion
    simoptions = setfieldifabsent(simoptions, 'maxAllowedxT', inf);
    
    simoptions = setfieldifabsent(simoptions, 'buoy', []);
    
    % Set up buoy simulation
    simoptions = buoysimsetup(simoptions.buoy, simoptions);
    
    % add zero WEC friction if not supplied
    design = setfieldifabsent(design, 'WECFriction', 0);

end