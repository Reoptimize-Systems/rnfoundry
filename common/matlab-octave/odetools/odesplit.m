function [results, tspan] = odesplit(odefcn, odeevfcn, tspan, y0, odeoptions, spfcn, varargin)
% odesplit: splits the evaluation of a system of differential equations
% into manageable chunks
%
% Syntax
%
% [results] = odesplit(odefcn, odeevfcn, tspan, y0, odeoptions, spfcn)
% 
% [results, tspan] = odesplit(odefcn, odeevfcn, tspan, y0, odeoptions, spfcn)
%
% [...] = odesplit(odefcn, odeevfcn, tspan, y0, odeoptions, spfcn, 'Parameter', value)
%
% Input
%
%   odefcn - this is a function handle to the ode solver function to be
%            used to find the solution, it must have the same calling
%            syntax as the built-in matlab solvers, e.g. ode45, ode15s etc.
% 
%   odeevfcn - this is the ode evaluation function, it must be of the form:
%
%                        dydt = odefun(t,y, varargin)
%
%              Additional input arguments are passed in by the
%              parameter-value pair 'OdeArgs' (see below).
% 
%   tspan - (1x2) Vector specifying the interval of integration. The solver
%           imposes the initial conditions at tspan(1), and integrates from 
%           tspan(1) to tspan(2), splitting the integration at the points
%           specified in the parameter-value pairs denoted by 'Blocks' and
%           'MaxDuration'(see below).
% 
%   y0 - Vector of initial conditions for the problem 
% 
%   odeoptions - Structure of optional parameters that change the default
%                integration properties. See matlab ode documentation and
%                function odeset.
% 
%   spfcn - Split point function. A function handle to the function to be 
%           called at each break in the integration. This function must
%           be able to be called with following syntax:
%
%           results = spfcn(flag, results, sol, arg0, arg1, ... , argn)
%
%           Where 'flag' is a scalar value indicating at what point in the
%           process of evaluation odesplit is at. The following values of
%           'flag' are possible:
%           
%           0 - this indicates odesplit is in the initialisation stageat
%               this point spfcn should return an initialized variable
%               (results) which will subsequently be passed in to it at
%               each split point in the integration along with the ode
%               solution structure and any additional arguments specified
%               in the p-v pair denoted by 'spfcnArgs' (See below).
%
%           1 - this indicates odesplit is calling the spfcn to perform any
%               processing required after a block of time has been
%               assessed.
%
%           2 - this indicates a terminal event has occured and odesplit is
%               calling the spfcn with the results of the block of time up
%               until this terminal event. 
% 
%           The state of any desired outputs from the ode solution can then
%           be maintained in the results variable for the duration of the
%           integration process. The sol structure is the normal ode
%           solution structure returned by the matlab solvers
%
%           sol.x: Steps chosen by the solver.
% 
%           sol.y: Each column sol.y(:,i) contains the solution at
%           sol.x(i).
% 
%           sol.solver: Solver name.
% 
%           If you specify the Events option and events are detected, sol
%           also includes these fields:
% 
%           sol.xe: Points at which events, if any, occurred. sol.xe(end)
%           contains the exact point of a terminal event, if any.
% 
%           sol.ye: Solutions that correspond to events in sol.xe.
% 
%           sol.ie: Indices into the vector returned by the function
%           specified in the Events option. The values indicate which event
%           the solver detected.
%
% odesplit will also detect terminal events and end the integration when
% these occur.
%
% odesplit also accepts the following optional parameter/value pairs (and
% therefore requires John D'Erico's parse_pv_pairs function):
%
%   'Blocks', scalar - the number of sections into which the time span
%   specified in 'tspan' should be split for evaluation. 
%
%   'MaxDuration', scalar - if a maximum duration is specified (in
%   seconds) then the ode is solved for blocks of time of length 'duration'
%   until tspan(2) is reached. The final section of time in this case may
%   be shorter than 'duration'. If you specify a MaxDuration longer than
%   tspan, the ode system will be solved in a single block.
%
%   If neither 'Blocks' nor 'MaxDuration' is specified then the time in
%   tspan is split into 2 equal sections. If both 'Blocks' and 'Maxduration'
%   are specified for some reason, MaxDuration takes precedence.
%
%   'OdeArgs', cell array - This is a cell array of additional arguments to
%   be passed to the function in 'odeevfcn'. 'odeevfcn' will be called by
%   the solver function in the following way:
%
%   @(t, y) feval(odeevfcn, t, y, OdeArgs{:})
%
%   'spfcnArgs', cell array - This is a cell array of additional arguments
%   to be passed to spfcn in addition to the results and solution
%   structures.
%
%   'ManageMemory', boolean - if true, odesplit will attempt to restart the
%   solution of the ode should it run out of memory, but using a greater
%   number of sections. By default odesplit will take two attempts
%   (including the initial attempt) and increase the number of sections by
%   a factor of 20 on the second run. I.e.if the number of blocks was
%   initially chosen to be 10, 200 blocks will be used for the second
%   attempt. Both of these aspects can also be controlled with the
%   'MaxAttempts' and 'BlockMultiplier' options. If the solution fails for
%   on the final attempt due to insufficient memory, the error is rethrown.
%   If the solution falis for any reason other than insufficient memory,
%   the error is immediately rethrown. Memory management is achieved by
%   wrapping the solution code in a try,catch block. If ManageMemory is set
%   to false the code is not solved in the Try,catch block. The default
%   value for ManageMemory is false.
%
%   'MaxAttempts', scalar - determines the number of times the solution
%   should be attempted in case of out of memory errors. The value of
%   MaxAttempts includes the initial attempt. On each attempt, the number
%   of sections into which the solution period is split is multiplied by a
%   factor. The default value of this factor is 20, but can be set by the
%   parameter 'BlockMultiplier'.
%
%   'BlockMultiplier', scalar - The multiplier for the number of sections
%   into which the solution period is split at each iteration in
%   'MaxAttempts'.
%
%   'Verbose', boolean - if true some information about the progress of the
%   solution is displayed, particularly regarding memory management
%   iterations if that option has been selected.
%
% Output
%
%   results - The final value of the variable/structure returned by spfcn
%             after the final solver run.
%
%   tspan - (1 x 2) vector, this contains the actual time span of the
%           simulation. Typically this will be the same as passed in, but
%           tspan(2) may be reduced if a terminal event occured during
%           simulation.
%
% Example
%
% Create the ode function described in Example 1 of the documentation for
% ode45 (doc ode45)
%
%   function dy = rigid(t,y)
%       dy = zeros(3,1);    % a column vector
%       dy(1) = y(2) * y(3);
%       dy(2) = -y(1) * y(3);
%       dy(3) = -0.51 * y(1) * y(2);
%
% Create a function to be called at each split point
% 
%    function results = rigidspfcn(flag, results, sol)
%        if flag == 0
%            results.iters = 0;
%        else
%            save(['rigidodesave_', int2str(results.iters), '.mat'], 'sol');
%            results.iters = results.iters + 1;
%        end
% 
% Set the options etc. and call odesplit, splitting the computation into
% three blocks
%
%   odeoptions = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5]);
%
%   odefcn = @ode45;
%   odeevfcn = @rigid;
%   tspan = [0, 30];
%   y0 = [0 1 1];
%   spfcn = @rigidspfcn;
%
%   [results, tspan] = odesplit(odefcn, odeevfcn, tspan, ...
%                           y0, odeoptions, spfcn, 'Blocks', 3)
%
%
%   See also 
%     ODE solvers:         ODE23, ODE45, ODE113, ODE15S, ODE23S, ODE23T, ODE23TB
%     options handling:    ODESET, ODEGET
%     evaluating solution: DEVAL
%     function handles:    FUNCTION_HANDLE

% Author: Richard Crozier
% Created: 09-12-2010 (That's december you silly Americans by the way!)
%

    Inputs.MaxDuration = 0;
    Inputs.Blocks = 2;
    Inputs.OdeArgs = {};
    Inputs.SplitPointFcnArgs = {};
    Inputs.ManageMemory = false;
    Inputs.BlockMultiplier = 20;
    Inputs.MaxAttempts = 2;
    Inputs.Verbose = false;
    
    Inputs = parse_pv_pairs(Inputs, varargin);
    
    % If a string is passed in as the evaluation function, convert this to
    % a function handle
    if ischar(odeevfcn)
        odeevfcn = str2func(odeevfcn);
    end
    
    if Inputs.MaxDuration == 0
        duration = (tspan(2)-tspan(1)) / Inputs.Blocks;
    elseif Inputs.MaxDuration < 0
        error('ODE block duration cannot be negative');
    else
        if Inputs.MaxDuration > (tspan(2) - tspan(1))
            Inputs.Blocks = 1;
        else
            Inputs.Blocks = ceil((tspan(2) - tspan(1)) / Inputs.MaxDuration);
        end
        
        duration =  (tspan(2)-tspan(1)) / Inputs.Blocks;
    end
    
    % first determine the new time intervals
    blocktspan(1) = tspan(1);
    blocktspan(2) = blocktspan(1) + duration;
    
    % Initialize the split point output function by calling it with the
    % initialisation flag
    flag = 0;
    results = feval(spfcn, flag, [], [], Inputs.OdeArgs{:});
    
    % The user may have specified the 'Refine' option for the ode solver, we
    % will use this value to ensure a sensible initial step when restarting
    % the solver on each loop 
    if isempty(odeoptions.Refine)
        refine = 1;
    else
        refine = odeoptions.Refine;
    end
    
    attempts = 1;
    blockcount = 0;
    
    if Inputs.Verbose
        % reset textprogressbar
        textprogressbar([]);
        % initialise textprogressbar
        textprogressbar('Split Ode Evaluation Progress: ');
        textprogressbar(100 * blockcount/Inputs.Blocks);
    end
     
    % for each block of time solve the ode by calling the appropriate
    % function, calling the split point function at the end of each 
    % block
    while blocktspan(2) <= tspan(2)
        
        istermevent = false;
        
        if Inputs.ManageMemory
            % If the ManageMemory value is true we repeat the entire ode
            % solution but with an increasing number of blocks if we run
            % out of memory until the max number of attempts is reached.
            try
                
                sol = odefcn(odeevfcn, blocktspan, y0, odeoptions, Inputs.OdeArgs{:});
                
                blockcount = blockcount + 1;
                
                if Inputs.Verbose
                    textprogressbar(100 * blockcount/Inputs.Blocks);
                end
                
            catch
                % for compatibility with octave, get the error by calling
                % lasterror
                ME = lasterror;

                % if the error was as a result of being out of memory and
                % we aren't out of attempts, then split up into smaller
                % blocks and try again.
                if ~isempty(strfind(ME.identifier, 'nomem')) && attempts < Inputs.MaxAttempts
                    
                    clear sol;
                    attempts = attempts + 1;
                    Inputs.Blocks = ceil(Inputs.Blocks * Inputs.BlockMultiplier);
                    duration =  (tspan(2)-tspan(1)) / Inputs.Blocks;
                    blocktspan(1) = tspan(1);
                    blocktspan(2) = blocktspan(1) + duration;
                    flag = 0;
                    results = feval(spfcn, flag, [], [], Inputs.OdeArgs{:});
                    blockcount = 0;
                    
                    if Inputs.Verbose
                        textprogressbar('\nEvaluation terminated due to lack of memory.');
                        fprintf(1, 'odesplit: Ran out of memory, increasing number of blocks to %d and beginning attempt number %d\n', Inputs.Blocks, attempts);
                        textprogressbar('Split Ode Evaluation Progress: ');
                        textprogressbar(100 * blockcount/Inputs.Blocks);
                    end
                    
                    continue;
                else
                    % An unexpected error happened, or we have reached the
                    % limit of allowed attempts and still ran out of memory
                    if Inputs.Verbose
                        textprogressbar('\nEvaluation terminated due to non memory related error.');
                        fprintf(1, 'odesplit: Ran out of memory using %d blocks on attempt number %d\n', Inputs.Blocks, attempts);
                    end
                    rethrow(ME);
                end

            end

        else
            sol = odefcn(odeevfcn, blocktspan, y0, odeoptions, Inputs.OdeArgs{:});
            
            blockcount = blockcount + 1;

            if Inputs.Verbose
                textprogressbar(100 * blockcount/Inputs.Blocks);
            end
        end
        
        % Now check for terminal events being reported in the solution
        % structure
        if isfield(sol, 'ie') && ~isempty(sol.ie)
            for ind = 1:numel(sol.ie)
                % odesplit expects the error function to take the same
                % additional arguments as the evaluation function
                [value,isterminal] = feval(odeoptions.Events, sol.xe(ind), sol.ye(:, ind), Inputs.OdeArgs{:});

                if isterminal(sol.ie(ind))
                    istermevent = true;
                    % set the flag to 2 to warn the split results
                    % evaluation function that a terminal event has been
                    % detected
                    flag = 2;
                    results = feval(spfcn, flag, results, sol, Inputs.SplitPointFcnArgs{:});
                end
            end
        else
           flag = 1;
           % Call the split point function with the solution results
           results = feval(spfcn, flag, results, sol, Inputs.SplitPointFcnArgs{:});
        end
        
        if sol.x(end) >= tspan(2) 
            
            % return as we're finished   
            if Inputs.Verbose
                textprogressbar('\nEvaluation complete.');
            end
            return;
            
        elseif istermevent
            
            % return as we hit a terminal event in the ode eval 
            if Inputs.Verbose
                textprogressbar('\nTerminal ODE event detected, evaluation halted.');
            end
            return;
            
        else
            
            % set the initial conditions of the system to be the final
            % condition of the previous block
            y0 = sol.y(:,end);

            % Note that starting from the end point of the last solution
            % will duplicate this time point in the output solution
            % received at the next iteration, i.e. if sol.x(end) = 10, at
            % the end of this iteration, the next iteration will also
            % include this time point so that sol.x(1) = 10. This may need
            % to be accounted for in your spfcn.
            blocktspan(1) = sol.x(end);
            
            % set the final time of the next block
            blocktspan(2) = min(blocktspan(1) + duration, tspan(2));
            
            % A good guess of a valid first time step is the length of
            % the last valid time step, so use it for faster computation.
            odeoptions = odeset(odeoptions, 'InitialStep', (sol.x(end) - sol.x(end-refine))/10);

        end
        
        % clear the solution to free up memory for the next calculation
        clear sol;
        
    end
    

    
end
    
    
