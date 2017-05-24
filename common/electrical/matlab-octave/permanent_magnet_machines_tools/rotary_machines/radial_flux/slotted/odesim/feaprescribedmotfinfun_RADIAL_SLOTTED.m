function [design, simoptions] = feaprescribedmotfinfun_RADIAL_SLOTTED (design, simoptions)
% set up an ode simulation of a slotted radial flux electrical machine
% based on full finite element analysis at every time step
%
% Syntax
%
% [design, simoptions] = feaprescribedmotfinfun_RADIAL_SLOTTED (design, simoptions)
%
%

    [design, simoptions] = prescribedmotfinfun_ROTARY (design, simoptions, @finfun_RADIAL_SLOTTED);

    [flux_linkage, ~, ~] = feaode_RADIAL_SLOTTED (design, simoptions, 0, simoptions.ODESim.SolutionComponents.PhaseCurrents.InitialConditions(:));

    simoptions.ODESim.SolutionComponents = setfieldifabsent (simoptions.ODESim.SolutionComponents, ...
                                          'PhaseFluxLinkages', ...
                                          struct ('InitialConditions', flux_linkage(:)', ...
                                                  'AbsTol', repmat (0.05*max(flux_linkage), 1, design.Phases)));
    simoptions.ODESim.Vectorized = 'on';
    
    % complete the setup with the common function
    [design, simoptions] = feaprescribedmotfinfun_AM (design, simoptions);
    
end

