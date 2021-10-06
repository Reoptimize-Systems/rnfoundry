function [design, simoptions] = feaprescribedmotfinfun_TM_SLOTLESS (design, simoptions)
% set up an ode simulation of a linear slotless tubular electrical machine
% based on full finite element analysis at every time step
%
% Syntax
%
% [design, simoptions] = feaprescribedmotfinfun_TM_SLOTLESS (design, simoptions)
%
%

    [design, simoptions] = prescribedmotfinfun_linear (design, simoptions, @finfun_TM_SLOTLESS);

    [flux_linkage, ~, ~] = feaode_TM_SLOTLESS (design, simoptions, ...
        0, simoptions.ODESim.SolutionComponents.PhaseCurrents.InitialConditions(:));

% [flux_linkage, ~, ~] = feaode_TM_SLOTLESS (design, simoptions, ...
%     0, zeros (design.Phases, 1));

% [flux_linkage, ~, ~, ~, ~] = ...
%         slmflmachineodesim_AM (design, simoptions, ...
%                                0, 0, ...
%                                1.0, 0, ...
%                                simoptions.ODESim.SolutionComponents.PhaseCurrents.InitialConditions(:) );

    simoptions.ODESim.SolutionComponents = setfieldifabsent (simoptions.ODESim.SolutionComponents, ...
                                          'PhaseFluxLinkages', ...
                                          struct ('InitialConditions', flux_linkage(:)', ...
                                                  'AbsTol', repmat (0.05*max(flux_linkage), 1, design.Phases)));
    simoptions.ODESim.Vectorized = 'on';
    
    % complete the setup with the common function
    [design, simoptions] = feaprescribedmotfinfun_AM (design, simoptions);
    
end

