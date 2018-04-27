function simoptions = simsetup_linear (design, PreProcFcn, PostPreProcFcn, varargin)
% TODO: help for simsetup_linear

    Inputs.TSpan = [];
    Inputs.Velocity = [];
    Inputs.TVelocity = [];
    Inputs.PoleCount = [];
    Inputs.RampPoles = [];
    Inputs.EvalFcn = 'prescribedmotodeforcefcn_linear'; 
    Inputs.PostSimFcn = 'prescribedmotresfun_linear'; 
    Inputs.ForceFcn = 'forcefcn_linear_pscbmot';
    Inputs.MinPointsPerPole = 10;
    Inputs.simoptions = struct();
    Inputs.ForceFcnArgs = {};
    
    Inputs = parse_pv_pairs (Inputs, varargin);
    
    simoptions = Inputs.simoptions;
    
    simoptions = setfieldifabsent (simoptions, 'ODESim', struct ());
    
    % strip existing sim spec if present
    simoptions = rmiffield (simoptions, 'xT');
    simoptions = rmiffield (simoptions, 'vT');
    simoptions = rmiffield (simoptions, 'drivetimes');

    if ~isempty (Inputs.Velocity)

        if isempty (Inputs.TSpan)
            if ~isempty (Inputs.PoleCount)
                Inputs.TSpan = [0, Inputs.PoleCount * design.PoleWidth / Inputs.Velocity];
            else
                error('SIMSETUP_ROTARY:notspan', ...
                    'If supplying constant Velocity  and not a pole corss count you must also supply the time span of the simulation.');
            end
        end

        ninterppoints = 10;

        simoptions.drivetimes = linspace (Inputs.TSpan(1), Inputs.TSpan(2), ninterppoints);

        simoptions.vT = repmat (Inputs.Velocity, 1, ninterppoints);

    elseif ~isempty (Inputs.TVelocity)

        if ~isfield (simoptions, 'drivetimes')

            simoptions.drivetimes = Inputs.TVelocity(:,1)';

        end

        if ~(isvector (simoptions.drivetimes) && all (diff (simoptions.drivetimes)>0))

            error('SIMSETUP_ROTARY:baddrivetimes', ...
                'The times supplied were not a monatonically increasing vector.');

        end

        simoptions.xt = Inputs.TVelocity(:,2)';

        simoptions.vt = Inputs.TVelocity(:,3)';

    else
        error('Insufficient information provided to set up simulation');
    end
    
    if isfield (simoptions, 'vT') && ~isfield (simoptions, 'xT')
        % v = dx / dt, so integrate to get the position
        simoptions.xT = cumtrapz(simoptions.drivetimes, simoptions.vT);
    end

    if ~isfield (simoptions.ODESim, 'EvalFcn') || isempty (simoptions.ODESim.EvalFcn)
        simoptions.ODESim.EvalFcn = Inputs.EvalFcn;
    end
    if ~isfield (simoptions.ODESim, 'PostSimFcn') || isempty (simoptions.ODESim.PostSimFcn)
        simoptions.ODESim.PostSimFcn = Inputs.PostSimFcn;
    end
    if ~isfield (simoptions.ODESim, 'PreProcFcn') || isempty (simoptions.ODESim.PreProcFcn)
        simoptions.ODESim.PreProcFcn = PreProcFcn;
    end
    if ~isfield (simoptions.ODESim, 'PostPreProcFcn') || isempty (simoptions.ODESim.PostPreProcFcn)
        simoptions.ODESim.PostPreProcFcn = PostPreProcFcn;
    end
    if ~isfield (simoptions.ODESim, 'ForceFcn') || isempty (simoptions.ODESim.ForceFcn)
        simoptions.ODESim.ForceFcn = Inputs.ForceFcn;
    end
    if ~isfield (simoptions.ODESim, 'ForceFcnArgs') || isempty (simoptions.ODESim.ForceFcnArgs)
        simoptions.ODESim.ForceFcnArgs = Inputs.ForceFcnArgs;
    end
    
    % if an initial ramp up in speed has been specified, construct it
    if ~isempty (Inputs.RampPoles) && Inputs.RampPoles > 0
        
        % add a linear speed ramp up over the specified number of poles,
        % typically to reduce the starting currents due to inductance
        nramppoles = Inputs.RampPoles;
        rampa = simoptions.vT(1)^2 / (2 * nramppoles * design.PoleWidth);
        rampTmax = simoptions.vT(1) / rampa;
        rampT = linspace (0, rampTmax, 15);
        rampvT = rampa .* rampT;
        rampxT = 0.5 * rampa .* rampT.^2;

        simoptions.vT = [rampvT(1:end-1), simoptions.vT];
        simoptions.xT = [rampxT(1:end-1), simoptions.xT + rampxT(end)];
        simoptions.drivetimes = [rampT(1:end-1), simoptions.drivetimes + rampT(end)];

        simoptions.ODESim.TimeSpan = simoptions.drivetimes([1, end]);

        simoptions.ODESim = setfieldifabsent (simoptions.ODESim, 'MaxStep', ...
            (simoptions.ODESim.TimeSpan(end) - simoptions.ODESim.TimeSpan(end-1)) / (Inputs.MinPointsPerPole * Inputs.RampPoles) );
        
        simoptions.ODESim.InitialStep = simoptions.ODESim.MaxStep;
        
    end
    
    simoptions.ODESim.TimeSpan = [simoptions.drivetimes(1), simoptions.drivetimes(end)];
    
    simoptions.ODESim = setfieldifabsent (simoptions.ODESim, 'MaxStep', maxstep (design, simoptions, Inputs.MinPointsPerPole));
    
    % construct a piecewise polynomial interpolation of the position
    % and velocity data
    simoptions.pp_xT = interp1 (simoptions.drivetimes, simoptions.xT, 'pchip', 'pp');
    simoptions.pp_vT = interp1 (simoptions.drivetimes, simoptions.vT, 'pchip', 'pp');
    
    % we will make the minimum phase current of interest that which
    % generates a power of 10W per coil at 1m/s, or a current density of
    % 0.1 A/mm^2 in the winding, whichever is less
    minIofinterest = min(design.ConductorArea * 0.1e6, ...
                         (10 / (design.Maxdlambdadx)) ) * design.Branches;
                     
    switch Inputs.EvalFcn
        
        
        case { 'prescribedmotode_linear', ...
               'prescribedmotodeforcefcn_linear', ...
               'feaprescribedmotodeforcefcn_linear' }
            
            % create the phase current solution component specification
            simoptions.ODESim.SolutionComponents = setfieldifabsent (simoptions.ODESim.SolutionComponents, ...
                                                  'PhaseCurrents', ...
                                                  struct ('InitialConditions', zeros (1, design.Phases), ...
                                                          'AbsTol', repmat (minIofinterest, 1, design.Phases) ) ...
                                                                    );
                                                                
%         case { 'prescribedmotodetorquefcn_dqactiverect_linear' }
% 
%             % create the phase current solution component specification
%             simoptions.ODESim.SolutionComponents = setfieldifabsent (simoptions.ODESim.SolutionComponents, ...
%                                                   'DQPhaseCurrents', ...
%                                                   struct ('InitialConditions', zeros (1, 2), ...
%                                                           'AbsTol', repmat (minIofinterest, 1, 2) ) ...
%                                                                     );
                                                                
        otherwise
            
            
            
    end
    
    
end

function dT = maxstep (design, simoptions, pperpole)
% choose a suitible max allowed time step

    maxvT = max(simoptions.vT);
    
    dT = design.PoleWidth / maxvT / pperpole;

end

