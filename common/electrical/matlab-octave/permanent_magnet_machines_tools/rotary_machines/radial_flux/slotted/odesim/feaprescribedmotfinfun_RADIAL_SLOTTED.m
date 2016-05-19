function [design, simoptions] = feaprescribedmotfinfun_RADIAL_SLOTTED(design, simoptions)

    [design, simoptions] = prescribedmotfinfun_ROTARY(design, simoptions, @finfun_RADIAL_SLOTTED);

    [flux_linkage, ~, ~] = feaode_RADIAL_SLOTTED (design, simoptions, 0, simoptions.ODESim.SolutionComponents.PhaseCurrents.InitialConditions);

    simoptions.ODESim.SolutionComponents = setfieldifabsent (simoptions.ODESim.SolutionComponents, ...
                                          'PhaseFluxLinkages', ...
                                          struct ('InitialConditions', flux_linkage(:)' ));
                               
	% set up the numerical derivative calculation of the flux linkage using
	% an odederiv object
    flodederiv = odederiv (0, simoptions.ODESim.SolutionComponents.PhaseFluxLinkages.InitialConditions');
    
	simoptions.ODESim.SolutionComponents.PhaseFluxLinkages.NumericalDerivatives = ...
        flodederiv;
    
    % we also create an output function for these components to update the
    % derivatives after each time step
    simoptions.ODESim.SolutionComponents.PhaseFluxLinkages.OutputFcn = ...
        @flodederiv.outputfcn;
    
    simoptions.ODESim = setfieldifabsent (simoptions.ODESim, 'OutputFcn', 'odesimoutputfcns_AM');
    
end