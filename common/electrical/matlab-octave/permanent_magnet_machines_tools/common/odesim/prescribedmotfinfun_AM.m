function [design, simoptions] = prescribedmotfinfun_AM (design, simoptions, finfun)

    % ensure any existing ode absolute tolerances are stripped so only the
    % values set in finfun, if any, are used
    simoptions = rmiffield(simoptions, 'abstol');

    if ~all(isfield(design, {'slm_fluxlinkage'}))
        % In this case we assume we have not already run the finalisation
        % code on this design and must do so
        [design, simoptions] = feval(finfun, design, simoptions);
    end

    % create the phase current solution component specification
    simoptions.ODESim.SolutionComponents = setfieldifabsent (simoptions.ODESim.SolutionComponents, ...
                                          'PhaseCurrents', ...
                                          struct ('InitialConditions', zeros (1, design.Phases) ));

end