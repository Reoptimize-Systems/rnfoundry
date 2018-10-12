function simoptions = simsetup_AM (design, simoptions, evalfcn, varargin)
% performs simulation setup operations common to all mahchine tyels

    Inputs.ForceAddPhaseCurrentODESolutionComps = false;
    Inputs.NoAddPhaseCurrentODESolutionComps = false;
    Inputs.ForceAddDQCurrentODESolutionComps = false;
    Inputs.NoAddDQCurrentODESolutionComps = false;
    Inputs.PhaseCurrentSolCompFcns = {'dummy'};
    Inputs.DQCurrentSolCompFcns = {'dummy'};
    
    Inputs = parse_pv_pairs (Inputs, varargin);
    
    check.isLogicalScalar (Inputs.ForceAddPhaseCurrentODESolutionComps, 'ForceAddPhaseCurrentODESolutionComps', true);
    check.isLogicalScalar (Inputs.NoAddPhaseCurrentODESolutionComps, 'NoAddPhaseCurrentODESolutionComps', true);
    check.isLogicalScalar (Inputs.ForceAddDQCurrentODESolutionComps, 'ForceAddDQCurrentODESolutionComps', true);
    check.isLogicalScalar (Inputs.NoAddDQCurrentODESolutionComps, 'NoAddDQCurrentODESolutionComps', true);
    
    assert (~(ForceAddPhaseCurrentODESolutionComps && NoAddPhaseCurrentODESolutionComps), ...
        'Both ForceAddPhaseCurrentODESolutionComps and NoAddPhaseCurrentODESolutionComps are true.' );
    assert (~(ForceAddDQCurrentODESolutionComps && NoAddDQCurrentODESolutionComps), ...
        'Both ForceAddDQCurrentODESolutionComps and NoAddDQCurrentODESolutionComps are true.' );
    assert (iscellstr (Inputs.PhaseCurrentSolCompFcns), ...
        'PhaseCurrentSolCompFcns must be a cell array of character vectors' );
    assert (iscellstr (Inputs.PhaseCurrentSolCompFcns), ...
        'DQCurrentSolCompFcns must be a cell array of character vectors' );
    
    % we will make the minimum phase current of interest that which
    % generates a power of 10W per coil at 1m/s, or a current density of
    % 0.1 A/mm^2 in the winding, whichever is less
    if isfield (design, 'ConductorArea')
        if isfield (design, 'Maxdlambdadx')
            minIofinterest = min ( design.ConductorArea * 0.1e6, ...
                                   (10 / (design.Maxdlambdadx)) ) ...
                                * design.Branches;
        else
            minIofinterest = design.ConductorArea * 0.1e6;
        end
    else
        minIofinterest = [];
    end
             
    if Inputs.ForceAddPhaseCurrentODESolutionComps
        % add the ODE solution components specification for full multiphase
        % current simulation
        
        simoptions.ODESim = setfieldifabsent (simoptions.ODESim, 'SolutionComponents', struct ());

        % create the phase current solution component specification
        simoptions.ODESim.SolutionComponents = setfieldifabsent ( simoptions.ODESim.SolutionComponents, ...
                                              'PhaseCurrents', ...
                                              struct ('InitialConditions', zeros (1, design.Phases) ) ...
                                                                );

        if ~isempty (minIofinterest)
            simoptions.ODESim.SolutionComponents.PhaseCurrents.AbsTol = repmat (minIofinterest, 1, design.Phases);
        end
        
    elseif Inputs.ForceAddDQCurrentODESolutionComps
        % add the ODE solution components specification for dq current
        % simulation
        
        simoptions.ODESim = setfieldifabsent (simoptions.ODESim, 'SolutionComponents', struct ());

        % create the phase current solution component specification
        simoptions.ODESim.SolutionComponents = setfieldifabsent (simoptions.ODESim.SolutionComponents, ...
                                              'DQPhaseCurrents', ...
                                              struct ('InitialConditions', zeros (1, 2) ) ...
                                                                );

        if ~isempty (minIofinterest)
            simoptions.ODESim.SolutionComponents.PhaseCurrents.AbsTol = repmat (minIofinterest, 1, 2);
        end
        
    else
        
        switch evalfcn

            case Inputs.PhaseCurrentSolCompFcns
                
                if ~Inputs.NoAddPhaseCurrentODESolutionComps

                    simoptions.ODESim = setfieldifabsent (simoptions.ODESim, 'SolutionComponents', struct ());

                    % create the phase current solution component specification
                    simoptions.ODESim.SolutionComponents = setfieldifabsent ( simoptions.ODESim.SolutionComponents, ...
                                                          'PhaseCurrents', ...
                                                          struct ('InitialConditions', zeros (1, design.Phases) ) ...
                                                                            );

                    if ~isempty (minIofinterest)
                        simoptions.ODESim.SolutionComponents.PhaseCurrents.AbsTol = repmat (minIofinterest, 1, design.Phases);
                    end
                
                end

            case Inputs.DQCurrentSolCompFcns
                
                if ~Inputs.NoAddDQCurrentODESolutionComps

                    simoptions.ODESim = setfieldifabsent (simoptions.ODESim, 'SolutionComponents', struct ());

                    % create the phase current solution component specification
                    simoptions.ODESim.SolutionComponents = setfieldifabsent (simoptions.ODESim.SolutionComponents, ...
                                                          'DQPhaseCurrents', ...
                                                          struct ('InitialConditions', zeros (1, 2) ) ...
                                                                            );

                    if ~isempty (minIofinterest)
                        simoptions.ODESim.SolutionComponents.PhaseCurrents.AbsTol = repmat (minIofinterest, 1, 2);
                    end
                
                end

            otherwise
                % do nothing
        end
    
    end
    
end

