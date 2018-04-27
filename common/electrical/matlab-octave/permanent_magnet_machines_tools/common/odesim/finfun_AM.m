function [design, simoptions] = finfun_AM(design, simoptions)
% performs common post-processing of machine design and simoptions prior to
% a dynamic ode simulation

% Copyright Richard Crozer, The University of Edinburgh

    % set any missing design properties
    
    % check the wire and coil properties are all present
    design = checkcoilprops_AM (design);
    
    % calculate the circuit properties
    [design, simoptions] = circuitprops_AM(design, simoptions);
    
    % calculate the coil wire conductor area
    design = setfieldifabsent(design, 'ConductorArea', circlearea (design.Dc/2));

    % coefficient of friction for translator (for legacy code)
    design = setfieldifabsent(design, 'mu_fT', 0);

    % coefficient of friction for armature (for legacy code)
    design = setfieldifabsent(design, 'mu_fA', 0);

    % coefficient of friction for effector
    design = setfieldifabsent(design, 'mu_fE', 0);

    % coefficient of friction for reactor
    design = setfieldifabsent(design, 'mu_fR', 0);

    % if absent define the Field direction to be 1, implying the magnetic
    % field moves with the effector, and the coil direction is the reverse
    % of this
    design = setfieldifabsent(design, 'FieldDirection', 1);

    % Set the coil positions (specified as a fraction of pole width for the
    % given Phases
    design = setfieldifabsent(design, 'CoilPositions', coilpos(design.Phases));
    
    % count the number of coils per phase
    design = setfieldifabsent(design, 'NCoilsPerPhase', design.CoilsPerBranch * design.Branches);
    
    % Set optional drag force proportional to velocity to zero
    design = setfieldifabsent(design, 'klineardrag', 0);
    
    % set the number of stages in the machine
    design = setfieldifabsent(design, 'NStages', 1);
    
    % Set the number of machines connected mechanically, 1 by default
    simoptions = setfieldifabsent(simoptions, 'NoOfMachines', 1);
    
    % set a flag which determines whether the results are saved to disk in
    % a split ode simulation
    simoptions = setfieldifabsent(simoptions, 'ODESim', struct());
    
%     simoptions.ODESim = setfieldifabsent(simoptions.ODESim, 'SolutionComponents', struct());
    simoptions.ODESim = odesim_set_recurse (simoptions.ODESim, 'SolutionComponents', struct());
    
    % set the number of outputs to use in the calcualtion of ode results
%     simoptions.ODESim = setfieldifabsent(simoptions.ODESim, 'ResultsTSkip', 1);
    simoptions.ODESim = odesim_set_recurse (simoptions.ODESim, 'ResultsTSkip', 1);
    
%     simoptions.ODESim = setfieldifabsent(simoptions.ODESim, 'SaveSplitResults', false);
    simoptions.ODESim = odesim_set_recurse (simoptions.ODESim, 'SaveSplitResults', false);
    
    % by default do not store all the results when getting splitting ode
    % solution into chunks
%     simoptions.ODESim = setfieldifabsent(simoptions.ODESim, 'StoreFullResults', false);
    simoptions.ODESim = odesim_set_recurse (simoptions.ODESim, 'StoreFullResults', false);
    
    simoptions = setfieldifabsent(simoptions, 'GetVariableGapForce', true);
    
    simoptions = setfieldifabsent (simoptions, 'TempDependantResistivity', false);
    
    simoptions = setfieldifabsent (simoptions, 'FreqDependantResistance', false);
    
    simoptions = setfieldifabsent(simoptions, 'basescorefcn', 'costscore_AM');

    if ~isfield (design, 'slm_fluxlinkage')
        
        if max(design.psilookup(1,:)) - min(design.psilookup(1,:)) <= 1
            % psilookup is provided over one pole, replicate the data to cover
            % two Poles
            flpos = [design.psilookup(1,:), design.psilookup(1,end) + design.psilookup(1,2:end)];  
            fl = [ design.psilookup(2,:), fliplr(design.psilookup(2,1:end-1))];

        else
            flpos = design.psilookup(1,:);  
            fl = design.psilookup(2,:);
        end

        % fit a periodic slm to the flux linkage against the normalised
        % positions
        design.slm_fluxlinkage = slmengine(flpos, fl, ...
                'EndCon', 'periodic', ...
                'knots', max(50, min(20, ceil(numel(design.psilookup)/1.5))), ...
                'Plot', 'off');
            
    end
    
    % calculate the percentage total harmonic distortion in the voltage
    % waveform produced
	design.VoltagePercentTHD = emfthd_AM(design.slm_fluxlinkage);
        
% 	if isfield(design, 'coggingforces')
%         
%         design.slm_coggingforce = slmengine(flpos, fl, ...
%                 'EndCon', 'periodic', ...
%                 'knots', min(20, ceil(numel(design.psilookup)/1.5)), ...
%                 'Plot', 'off');
%         
%     end
    
    % determine the maximum rate of change in flux linkage with
    % displacement for later use
    design.Maxdlambdadx = slmpar(design.slm_fluxlinkage, 'maxslope') / design.PoleWidth;
    
    % store the rms value of the flux linkage
    design.FluxLinkageRms = rms(slmeval(linspace(0,2,1000), design.slm_fluxlinkage, 0, false));
    
    design.FluxLinkageCoilPeak = slmpar(design.slm_fluxlinkage, 'maxfun');
    
    design.FluxLinkagePhasePeak = design.FluxLinkageCoilPeak * design.CoilsPerBranch;

    % by default we will add the phase currents to the list of solution
    % components to be solved
    simoptions = setfieldifabsent (simoptions, 'AddPhaseCurrentsComponents', false);
    
    if simoptions.AddPhaseCurrentsComponents   
        
        % we will make the minimum phase current of interest that which
        % generates a power of 10W per coil at 1m/s, or a current density of
        % 0.1 A/mm^2 in the winding, whichever is less
        minIofinterest = min(design.ConductorArea * 0.1e6, ...
                         (10 / (design.Maxdlambdadx)) ) * design.Branches;
                     
        % create the phase current solution component specification
        simoptions.ODESim.SolutionComponents = setfieldifabsent (simoptions.ODESim.SolutionComponents, ...
                                              'PhaseCurrents', ...
                                              struct ('InitialConditions', zeros (1, design.Phases), ...
                                                      'AbsTol', repmat (minIofinterest, 1, design.Phases) ) ...
                                                                );
    end
    
    % set infinite allowed deflection factor if none is supplied
    simoptions = setfieldifabsent(simoptions, 'maxAllowedDeflectionFactor', inf);
    
    % run a function to display the machine design if it is supplied
    if isfield(simoptions, 'DisplayDesignFcn')

        feval(simoptions.DisplayDesignFcn, design, simoptions);

    end
    
    % if the data is supplied, fit a curve to variation in air gap closing
    % force with variation in air gap
    if isfield(design, 'gforce') && isfield(design, 'gvar')
        if numel (design.gvar) > 2
            order = 2;
        else
            order = 1;
        end
        design.p_gforce = polyfitn(design.gvar, design.gforce, order);
        
    end
    
    design.PostPreProcessingComplete = true;
    
end


function ODESim = odesim_set_recurse (ODESim, fieldname, defaultval)
% recursively sets a default option for ODESim simulations and any Nested
% ODESim simulations
%

    if isfield (ODESim, 'NestedSim')
        ODESim.NestedSim = odesim_set_recurse (ODESim.NestedSim, fieldname, defaultval);
    end
    
    ODESim = setfieldifabsent (ODESim, fieldname, defaultval);

end
