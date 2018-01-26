function [T, Y, results, design, simoptions] = simulatemachine_AM(design, simoptions, varargin)
% performs a simulation of an electrical machine design and system using
% ode solvers
%
% Syntax
%
% [T, Y, results, design, simoptions] = simulatemachine_AM (design, simoptions, 'Parameter', 'Value')
%
%
% Description
%
% simulatemachine_AM performs preprocessing, runs an ode simulation,
% optionally with memory management, and post-processes the results using
% functions supplied in 'simoptions'. 
%
% Inputs
%
%  design - Stucture contining the specification of the design to be
%    simulated.
%
%  simoptions - structure used to control how the system is evaluated.
%   The simoptions structure must contain a substructure named ODESim. 
%
%   The ODESim structure can contain the following fields:
%
%   PreProcFcn : optional function handle or string containing the function
%     which will be run to generate data prior to running the simulation.
%     simfun will be passed the design and simoptions structure, and can
%     also be supplied with additional arguments by using the
%     'ExtraPreProcFcnArgs' parameter-value pair. The extra arguments must be
%     placed in a cell array. It must return two arguments which will
%     overwrite the design and simoptions variables.
%
%     The supplied function must have the calling syntax
%
%     [design, simoptions] = thefunction(design, simoptions, simfunarg1, simfunarg2, ...)
%
%     Where simfunarg1, simfunarg2, ... if supplied are the elements of the
%     a cell array, passed in using the Parameter-value pair
%     ExtraPreProcFcnArgs, e.g.
%
%     simulatemachine_AM(design, simoptions, 'ExtraPreProcFcnArgs', {1, 'another', [1,2;3,4]})
%
%   PostPreProcFcn : optional function handle or string containing a
%     function which will be run after simfun. finfun will also be passed
%     the design and simoptions structure, and can also be supplied with
%     additional arguments by using the 'ExtraPostPreProcFcnArgs'
%     parameter-value pair. The extra arguments must be placed in a cell
%     array. It must return two arguments which will overwrite the design
%     and simoptions variables.
%
%     The supplied function must have the calling syntax
%
%     [design, simoptions] = thefunction(design, simoptions, finfunarg1, finfunarg2, ...)
%   
%     Where finfunarg1, finfunarg2, ... if supplied are the elements of the
%     a cell array, passed in using the Parameter-value pair
%     ExtraPostPreProcFcnArgs, e.g.
%
%     simulatemachine_AM(design, simoptions, 'ExtraPostPreProcFcnArgs', {1, 'another', [1,2;3,4]})
%
%   EvalFcn : function handle or string containing the function which will
%     be evaluated by the ode solver routines to solve the system of
%     equations. see the ode solvers (e.g. ode45, ode15s) for further
%     information on how to create a suitible function.
%
%   SolutionComponents : Structure containing information on the various
%     components of the system to be solved. This structure will contain
%     fields with names corresponding to the solution components. Each of
%     these fields is then also a structure with information about each of
%     these solution components. The information is contained in the
%     following fields:
%
%     SolutionIndices : these are the indices of the ode system
%       coresponding to this component.
%
%     InitialConditions : this is the initial conditions for variables
%       associated with this component.
%
%     AbsTol : the absolute tolerances of the variables associated with
%       this  component. If AbsTol is not supplied for every solution
%       component it will be ignored for all components.
%
%     The initial conditions for the ODE solution and possibly the absolute
%     tolerances on all variables are assembled from the information in
%     this structure.
%
%     In addition, the field OutputFcn may be present. This contains a
%     function handle or string representing a function to be run on the
%     variables associated with the solution component after each
%     successful time step.
%
%     SolutionComponents can also contain a field named 'NestedSim' for the
%     purposes of setting up a multi-rate solution. If present, it contains
%     solution component information intended for a nested, higher time
%     step rate simulation (implemented using the ode.odesolver class),
%     running between the steps of the 'outer' or top level simulation. The
%     NestedSim sturcture is the same format as the top level
%     SolutionComponents structure and contains the solution component
%     information for the variables in the lower level simulation. The
%     NestedSim structure can also contain another NestedSim structure and
%     any number of levels are supported.
%
%   PostAssemblyFcn : function handle or string containing a function which
%     will be run after the solution components have been read and fully
%     assembled.
%
%   PostSimFcn : function handle or string containing a function which will
%     be run after the simulation has been completed by the ode solver.
%     resfun must take the T and Y matrices, as generated by the ode
%     solver, and the design and simoptions arguments in that order.
%     PostSimFcn can also be supplied with additional arguments by using
%     the 'ExtraPostSimFcnArgs' parameter-value pair. The extra arguments
%     must be placed in a cell array. It must return two arguments, one of
%     which is a results variable containing results of interest to the
%     user, the other of which overwrites the design variable.
%
%     The supplied function must have the calling syntax
%
%     [results, design] = htefunction(T, Y, design, simoptions, odearg1, odearg2, ..., resfunarg1, resfunarg1, ...);
%
%     Where odearg1, odearg2, ..., resfunarg1, resfunarg1, ... if supplied are the elements of the
%     two cell arrays, passed in using the Parameter-value pairs
%     ExtraEvalFcnArgs and ExtraPostSimFcnArgs respectively, e.g.
%
%     simulatemachine_AM ( design, simoptions, ...
%                          'ExtraEvalFcnArgs', {1, true}, ...
%                          'ExtraPostSimFcnArgs', {2, false} )
%
%   Solver : function handle or string specifying the ode solver to use,
%     if not supplied, for Matlab the default is 'ode15s', and for Octave
%     'ode2r'.
%
%   splitode : (scalar integer) if this field is present in the structure
%     it indicates that the evaluation of the system of differential
%     equations should be split into manageable chunks, useful for
%     long-running simulations which use substantial memory. The value of
%     splitode is the desired initial number of chunks into which
%     evaluation will be split. If the system runs out of memory during
%     simulation, the number of blocks will be increased and simulation
%     reattempted this will be attempted at most 4 times. If splitode is
%     present, the field spfcn must also be supplied, documented below.
%
%   SplitPointFcn : (string|function handle) if splitode is provided this
%     field must also be present which should contain a string or function
%     handle. This function will be called at each break in the integration
%     and must have the following syntax:
%
%     results = spfcn (flag, results, sol, design, simoptions)
%   
%     For further information on creating a split point function, see the
%     help for 'odesplit.m'.
%
%   OutputFcn : Top level OutputFcn to be run after each successful time
%     step. If not supplied, the function odesimoutputfcns_AM is used. This
%     function looks for a SolutionComponents structure in
%     simoptions.ODESim and runs the output function (if present) for each
%     solution component with the appropriate inputs.
%
%   Events : 
%
%   Vectorized : true/false flag indicating if the ODE solution function is
%     vectorized
%
%   MaxStep : The maximum allowed time step size
%
%   InitialStep : A suggested initial time step size for the ode solution
%
% Additional arguments are provided via Parameter-Value pairs, most of
% which are related to the contents of the simoptions fileds and are
% already described previously (such as the options to pass additional
% arguments to the simulation functions). These Parameter-Value pairs are:
%
%   'ExtraPreProcFcnArgs' - cell array of a additional arguments to pass to
%     simfun
%
%   'ExtraPostPreProcFcnArgs' - cell array of a additional arguments to pass to
%     funfun
%
%   'ExtraPostSimFcnArgs' - cell array of a additional arguments to pass to
%     resfun
%
%   'ExtraEvalFcnArgs' - cell array of a additional arguments to pass to
%     odeevfun
%
% In addition, the following parameter-value pairs may be supplied:
%
%   'Verbose' - true or false flag determining whether to print output
%     describing the progress of the simulation. Default is true.
%
% Output
%
% T - output time vector as produced by ode solver functions, e.g. ode45
%
% Y - output solution vector as produced by ode solver functions, e.g.
%   ode45
%
% results - the results as produced by the supplied function in resfun
%
% design - the design structure which may have been modified by the
%   supplied functions
%
% simoptions - the simoptions structure which may have been modified by the
%   supplied functions
%

% Copyright Richard Crozer 2012-2016

    % Do some parsing of optional input arguments
    Inputs.ExtraPreProcFcnArgs = {};
    Inputs.ExtraPostPreProcFcnArgs = {};
    Inputs.ExtraPostAssemblyFcnArgs = {};
    Inputs.ExtraPostSimFcnArgs = {};
    Inputs.ExtraEvalFcnArgs = {};
    Inputs.Verbose = true;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    % get the simululation pre-processing/data gathering function
    simfun = [];
    if isfield(simoptions.ODESim, 'PreProcFcn') 
        if ischar(simoptions.ODESim.PreProcFcn)
            simfun = str2func(simoptions.ODESim.PreProcFcn);
        else
            simfun = simoptions.ODESim.PreProcFcn;
        end
    end
    
    if ~isempty(simfun)
        % run the data gathering function for the design
        [design, simoptions] = feval(simfun, design, simoptions, Inputs.ExtraPreProcFcnArgs{:});
    end
    
    finfun = [];
    if isfield(simoptions.ODESim, 'PostPreProcFcn')
        if ischar(simoptions.ODESim.PostPreProcFcn)
            finfun = str2func(simoptions.ODESim.PostPreProcFcn);
        else
            finfun = simoptions.ODESim.PostPreProcFcn;
        end
    end

    if ~isempty(finfun)
        % now complete any further required modification or post-processing
        % required before running the ode simulation
        [design, simoptions] = feval(finfun, design, simoptions, Inputs.ExtraPostPreProcFcnArgs{:});
    end

    % if solution components info is provided construct the initial
    % conditions and tolerances etc from them
    if isfield (simoptions.ODESim, 'SolutionComponents')
        simoptions = assemble_ode_components (simoptions);
    end
    
    PostAssemblyFcn = [];
    if isfield(simoptions.ODESim, 'PostAssemblyFcn')
        if ischar(simoptions.ODESim.PostAssemblyFcn)
            PostAssemblyFcn = str2func(simoptions.ODESim.PostAssemblyFcn);
        else
            PostAssemblyFcn = simoptions.ODESim.PostAssemblyFcn;
        end
    end
    
    if ~isempty(PostAssemblyFcn)
        % now complete any further required modification or post-processing
        % required before running the ode simulation
        [design, simoptions] = feval(PostAssemblyFcn, design, simoptions, Inputs.ExtraPostAssemblyFcnArgs{:});
    end
    
    odeargs = [{design, simoptions}, Inputs.ExtraEvalFcnArgs];
    
% ---- No modifications to simoptions below this point will be propogated to the ODE function

    % construct the various functions for the ode solvers. This involves
    % constructing anonymous functions to pass in the extra ode solver
    % arguments as appropriate
    [odeevfun, resfun, spfcn, eventfcns, outputfcn] = checkinputs_simulatemachine_AM (simoptions, odeargs);
    
    % finally simulate the machine using ode solver, first setting some
    % options
    if ~isfield(simoptions.ODESim, 'RelTol')
        simoptions.ODESim.RelTol = 2e-2;
    end

    odeoptions = odeset('RelTol', simoptions.ODESim.RelTol);
    
    % TODO: why? 
    % choose initial step size, this is done 
%     odeoptions = odeset(odeoptions, 'InitialStep', simoptions.ODESim.TimeSpan(end) / 1000);
    
    if isfield (simoptions.ODESim, 'AbsTol')
        odeoptions = odeset(odeoptions, 'AbsTol', simoptions.ODESim.AbsTol);
    end
    
    if isfield (simoptions.ODESim, 'InitialStep')
        odeoptions = odeset(odeoptions, 'InitialStep', simoptions.ODESim.InitialStep);
    end

    if isfield(simoptions, 'maxstep')
        odeoptions = odeset(odeoptions, 'MaxStep', simoptions.maxstep);
    elseif isfield (simoptions.ODESim, 'MaxStep')
        odeoptions = odeset(odeoptions, 'MaxStep', simoptions.ODESim.MaxStep);
    end

    if isfield(simoptions, 'events')
        odeoptions = odeset(odeoptions, 'Events', eventfcns);
    elseif isfield (simoptions.ODESim, 'Events')
        odeoptions = odeset(odeoptions, 'Events', eventfcns);
    end
    
    if isfield(simoptions.ODESim, 'OutputFcn')
        odeoptions = odeset(odeoptions, 'OutputFcn', outputfcn);
    end
    
    if isfield(simoptions.ODESim, 'Vectorized')
        odeoptions = odeset(odeoptions, 'Vectorized', simoptions.ODESim.Vectorized);
    end

    % select the solver to use, by default we choose stiff solvers a
    % generally the machine inductance circuit requires this.
    if isfield(simoptions.ODESim, 'Solver')
        if ischar(simoptions.ODESim.Solver)
            odefcn = str2func(simoptions.ODESim.Solver);
        elseif isa(simoptions.ODESim.Solver, 'function_handle')
            odefcn = simoptions.ODESim.Solver;
        end
    else
        if isoctave
            odefcn = @ode5r; %s @oders; %@ode23s; %@ode2r;
        else
            odefcn = @ode15s;
        end
    end
    
    if ischar(odeevfun)
        % convert evaluation function to function handle if it's a string
        odeevfun = str2func(odeevfun);
    end
    
    if isfield(simoptions.ODESim, 'Split') 
        
        if isfield(simoptions.ODESim, 'SplitPointFcn') && isa(spfcn, 'function_handle')

            % in this case we use the odesplit function to allow longer
            % simulations to be run, only extracting pertinent values
            [results, simoptions.ODESim.TimeSpan] = odesplit ( odefcn, ...
                odeevfun, ...
                simoptions.ODESim.TimeSpan, ...
                simoptions.ODESim.InitialConditions, ...
                odeoptions, ...
                spfcn, ...
                'SplitPointFcnArgs', {}, ...
                'OdeArgs', {}, ...
                'Blocks', simoptions.ODESim.Split, ...
                'BlockMultiplier', 10, ...
                'ManageMemory', true, ...
                'MaxAttempts', 4, ...
                'Verbose', true );

            % append the results to the Input arguments for the resfun
            Inputs.ExtraPostSimFcnArgs = [Inputs.ExtraPostSimFcnArgs, {results}];

            % set T and Y to empty matrices
            T = [];
            Y = [];
        else
            error('simoptions.ODESim contains Split field, but has missing or invalid SplitPointFcn handle');
        end
        
    else
        if Inputs.Verbose, fprintf(1, '\nBeginning ode solution\n'); end
        %tic
        %[T, Y] = odefcn(@(t, y) feval(odeevfun, t, y, design, simoptions, Inputs.odeargs{:}), simoptions.ODESim.TimeSpan, simoptions.ODESim.InitialConditions, odeoptions);
        [T,Y] = odefcn ( odeevfun, ...
                         simoptions.ODESim.TimeSpan, ...
                         simoptions.ODESim.InitialConditions, ...
                         odeoptions );
        %toc
        if Inputs.Verbose, fprintf(1, 'ode solution complete\n'); end
    end

    % Obtain any internally calculated results of interest
    [results, design] = feval ( resfun, ...
                                T, Y, design, simoptions, ...
                                Inputs.ExtraEvalFcnArgs{:}, ...
                                Inputs.ExtraPostSimFcnArgs{:} );
    
%     if isoctave
%         simoptions = func2str_simulatemachine_linear(simoptions);
%     end
    
end


function [odeevfun, resfun, spfcn, eventfcns, outputfcn] = checkinputs_simulatemachine_AM (simoptions, odeargs)

    odeevfun = [];
    if isfield(simoptions.ODESim, 'EvalFcn') 
        if ischar(simoptions.ODESim.EvalFcn)
            odeevfun = str2func(simoptions.ODESim.EvalFcn);
            odeevfun = @(t,y) odeevfun (t, y, odeargs{:});
        else
            odeevfun = simoptions.ODESim.EvalFcn;
        end
    end

    if isfield(simoptions.ODESim, 'PostSimFcn') 
        if ischar(simoptions.ODESim.PostSimFcn)
            resfun = str2func (simoptions.ODESim.PostSimFcn);
        else
            resfun = simoptions.ODESim.PostSimFcn;
        end
    end
    
    outputfcn = [];
    if isfield(simoptions.ODESim, 'OutputFcn') 
        if ischar(simoptions.ODESim.OutputFcn)
            outputfcn = str2func(simoptions.ODESim.OutputFcn);
            outputfcn = @(t,y,flag) outputfcn (t,y,flag, odeargs{:});
        else
            outputfcn = simoptions.ODESim.OutputFcn;
        end
    end
    
    spfcn = [];
    if isfield(simoptions.ODESim, 'SplitPointFcn') 
        if ischar(simoptions.ODESim.SplitPointFcn)
            spfcn = str2func(simoptions.ODESim.SplitPointFcn);
            spfcn = @(flag, results, sol) spfcn (flag, results, sol, odeargs{:});
        end
    end
    
    eventfcns = {};
    if isfield(simoptions, 'events')
        if iscell(simoptions.events)
            eventfcns = cell (size (simoptions.events));
            for i = 1:numel(simoptions.events)
                if ischar(simoptions.events{i})
                    eventfcns{i} = str2func(simoptions.events{i});
                    eventfcns{i} = @(t,y) eventfcns{i}(t,y,odeargs{:});
                else
                    eventfcns{i} = simoptions.events{i};
                end
            end
        else
            if ischar(simoptions.events)
                eventfcns = str2func(simoptions.events);
                eventfcns = @(t,y) eventfcns(t,y,odeargs{:});
            else
                eventfcns = simoptions.events;
            end
        end
    end
end

function simoptions = func2str_simulatemachine_linear(simoptions)

    if isfield(simoptions.ODESim, 'PreprocFcn') && isa(simoptions.ODESim.PreProcFcn, 'function_handle')
        simoptions.ODESim.PreProcFcn = func2str(simoptions.ODESim.PreProcFcn);
    end

    if isfield(simoptions.ODESim, 'PostPreProcFcn') && isa(simoptions.ODESim.PostPreProcFcn, 'function_handle')
        simoptions.ODESim.PostPreProcFcn = func2str(simoptions.ODESim.PostPreProcFcn);
    end

    if isfield(simoptions.ODESim, 'EvalFcn') && isa(simoptions.ODESim.EvalFcn, 'function_handle')
        simoptions.ODESim.EvalFcn = func2str(simoptions.ODESim.EvalFcn);
    end

    if isfield(simoptions.ODESim, 'PostSimFcn') && isa(simoptions.ODESim.PostSimFcn, 'function_handle')
        simoptions.ODESim.PostSimFcn = func2str(simoptions.ODESim.PostSimFcn);
    end
    
    if isfield(simoptions.ODESim, 'SplitPointFcn') && isa(simoptions.ODESim.SplitPointFcn, 'function_handle')
        simoptions.ODESim.SplitPointFcn = func2str(simoptions.ODESim.SplitPointFcn);
    end

    if isfield(simoptions, 'events')
        if iscell(simoptions.events)
            for i = 1:numel(simoptions.events)
                if isa(simoptions.events{i}, 'function_handle')
                    simoptions.events{i} = func2str(simoptions.events{i});
                end
            end
        else
            if isa(simoptions.events, 'function_handle')
                simoptions.events = func2str(simoptions.events);
            end
        end
    end
    
end

function simoptions = assemble_ode_components (simoptions)

    simoptions.ODESim = assemble_ode_comp_recurse (simoptions.ODESim);
    
    simoptions.ODESim = setfieldifabsent (simoptions.ODESim, 'OutputFcn', 'odesimoutputfcns_AM');

end

function odeoptstruct = assemble_ode_comp_recurse (odeoptstruct)
% function to recursively set up a set of nested ode simulations

    if isfield (odeoptstruct, 'NestedSim')
        % set up the simulation for the nested simulation
    	odeoptstruct.NestedSim = assemble_ode_comp_recurse (odeoptstruct.NestedSim);
    end

    compnames = fieldnames (odeoptstruct.SolutionComponents);

    odeoptstruct.InitialConditions = [];
    
    odeoptstruct.AbsTol = [];

    for ind = 1:numel (compnames)

        odeoptstruct.SolutionComponents.(compnames{ind}).SolutionIndices = ...
            (numel(odeoptstruct.InitialConditions) + 1) : ...
                numel(odeoptstruct.InitialConditions)+numel(odeoptstruct.SolutionComponents.(compnames{ind}).InitialConditions);

        odeoptstruct.InitialConditions = [ odeoptstruct.InitialConditions, ...
                odeoptstruct.SolutionComponents.(compnames{ind}).InitialConditions(:)' ];
            
        if isfield (odeoptstruct.SolutionComponents.(compnames{ind}), 'AbsTol')
            abstol = odeoptstruct.SolutionComponents.(compnames{ind}).AbsTol;
        else
            abstol = nan (size(odeoptstruct.SolutionComponents.(compnames{ind}).InitialConditions));
        end
        
        odeoptstruct.AbsTol = [ odeoptstruct.AbsTol, ...
                                     abstol(:)' ];

    end
    
    if any (isnan(odeoptstruct.AbsTol))
        odeoptstruct = rmfield (odeoptstruct, 'AbsTol');
        warning ('RENEWNET:simulatemachine_AM:badabstol', ...
                 ['AbsTol not supplied for all solution components. ALL AbsTol ', ...
                  'specifications have therefore been removed and will not be applied to ', ...
                  'the ode solution.'])
    end


end

