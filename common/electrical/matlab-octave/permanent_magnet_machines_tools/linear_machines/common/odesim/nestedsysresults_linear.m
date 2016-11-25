function results = nestedsysresults_linear (flag, results, sol, mc, evalfcn, design, simoptions)
% nestedsysresults_linear: accumulates data from the evaluation of a
% linear machine and heaving buoy simulation when using odesplit to split the
% evaluation time span into sections

    if flag == 0
        
        results.block = 1;
        
        % set up initial values for the electrical results
        results = splitodeelectricalresults_AM (design, simoptions, flag, ...
                            results, [], [], ...
                            simoptions.ODESim.NestedSim.SolutionComponents.PhaseCurrents.SolutionIndices, ...
                            design.NStages, ...
                            simoptions.ODESim.NestedSim.StoreFullResults );
        
        % Initialize the results structure
        results.minLongerPartLength = 0;
        
        results.interpdur = [];
        
        results.vRmax = 0;
        
        return;

    elseif flag == 1

        % Now obtain internally calculated values by recalling the function
        % with the results, first preallocating the arrays. This is necessary
        % as the ode solver used may take steps while choosing step sizes which
        % do not form part of the solution      
        odeinternals = oderesults (sol, mc, evalfcn);

        % get the maximum relative velocity
        results.vRmax = max(results.vRmax, max(abs(odeinternals.vE)));
        
        % Actual armature length that would be required to achieve the
        % output, this is not used in the cost calculations
        peakxT = max(abs(odeinternals.xE));

        % get the minum length the longer part must be to allow the range of
        % motion experienced by the translator
        results.minLongerPartLength = max(results.minLongerPartLength, 2 * peakxT + (design.PowerPoles * design.PoleWidth));

        % Determine some interesting design outputs
        results = splitodeelectricalresults_AM (design, simoptions, flag, results, sol, odeinternals, ...
                            simoptions.ODESim.NestedSim.SolutionComponents.PhaseCurrents.SolutionIndices, ...
                            design.NStages, ...
                            simoptions.ODESim.NestedSim.StoreFullResults );
        
                        
        if simoptions.ODESim.NestedSim.StoreFullResults
            if results.block == 1
                results.Time = sol.x';
                results.xE = odeinternals.xE;
                results.vE = odeinternals.vE;
            else
                results.Time = [results.Time; sol.x(2:end)'];
                results.xE = [ results.xE; odeinternals.xE(2:end) ];
                results.vE = [ results.vE; odeinternals.vE(2:end) ];
            end
        end
        
        results.block = results.block + 1;
        
%         % deal with any terminal events if these have occured
%         if isfield(sol, 'ie') && ~isempty(sol.ie)
%             
%             % append the terminal event info from the solution structure to
%             % the results
%             results.ie = sol.ie;
%             results.xe = sol.xe;
%             results.ye = sol.ye;
%             
%             % display a message dependent on the type of terminal event
%             % detected
%             fprintf(1, '\nTerminal event occured at time t=%f: ', results.xe(1));
%             
%             for ind = 1:numel(results.ie)
%                 
%                 switch results.ie(ind)
% 
%                     case 1
%                         fprintf(1, 'Buoy moving too fast in heave, speed = %f\n', results.ye(2,ind));
%                     case 2
%                         fprintf(1, 'Buoy moving too fast in surge, speed = %f\n', results.ye(4,ind));
%                     case num2cell(3:(2+design.Phases))
%                         fprintf(1, 'Too much current in coil %d, current = %f\n', results.ie(ind)-2, results.ye(results.ie(ind)+2,ind));
%                     otherwise
%                         fprintf(1, 'Unknown termination reason, event solution vector ind: %d\n', results.ie(ind));
%                         
%                 end
%                 
%             end
% 
%         end
        
    end

end


function results = oderesults(sol, mc, odeevfcn, odeargs, skip, skipfields)
% oderesults: extracts results from an appropriately coded ode which were
% not variables of the integration.
%
% Syntax
%
% results = oderesults(T, Y, odeevfcn)
% results = oderesults(T, Y, odeevfcn, odeargs)
% results = oderesults(T, Y, odeevfcn, odeargs, skip)
%
% Description
%
% oderesults recalls a function designed to be used with the standard
% matlab ODE solvers (e.g. ode45, ode15s etc.) after a solution is complete
% in order to extract internally calculated values of interest generated
% during the simulation, but which are not returned by the solver (which
% generally only return the variables which are the subject of the
% integration). To achieve this with oderesults, the solution function must
% be coded appropriately to return certain information about the desired
% results when called in a certain way.
%
% sol - solution structure as produced by the ODE solver using the solution
%   function in 'odeevfcn'
%
% odeevfcn is an ode solution function that can be called with multiple
%   output arguments. Each output argument must return a scalar or vector
%   which represents a result of interest at each time in T.  It must also
%   be coded so that when called with no arguments it returns a cell array
%   of strings which will become the names of fields in a results
%   structure. The number of results of interest will be determined by the
%   number of strings in the cell array. If returning a vector, the length
%   of the vector on each call must remain the same as they will be
%   concatenated into a matrix.
%
% odeargs in an optional cell array of additional arguments to passed to the
%   ode solution function, cal be an empty cell array if no arguments are
%   desired, but the 'skip' input is to be used.
%
% skip is a scalar integer which determines how much of the solution to
%   recalculate. If skip is greater than 1, only every 'skip'th value is
%   recalculated.
%
% skipfields is an optional cell array of output field names to ignore when
%   populating the results structure
%
% See also: prallocresfcn.m
%

    if nargin < 4
        odeargs = {};
    end
    
    if nargin < 5
        skip = 1;
    end
    
    if nargin < 6
        skipfields = {};
    end

    % preallocate the results structure with fields containing arrays of
    % zeros
    [results, fn, nOut] = prallocresfcn (sol, mc, odeevfcn, odeargs, skip);
    
    if skip ~= 1
        % store the times the results values occur at if we are skipping
        % some in the results structure
        results.Tskip = sol.x(1:skip:length (sol.x))';
    end
    
    k = 0;
    
    outArgs = cell(1, nOut);
    
    % now recalculate the values and put them in the results structure
    for i = 1:skip:length(sol.x)

        k = k + 1;
        
        [outArgs{1:nOut}] = feval (odeevfcn, sol.x(i), sol.y(:,i), mc, 1, odeargs{:});
        
        for j = 1:nOut
            if ~ismember (fn{j}, skipfields)
                results.(fn{j})(k,:) = outArgs{j};
            end
        end
    
    end
    
end


function [results, fn, nOut] = prallocresfcn(sol, mc, odeevfcn, odeargs, skip)
% prallocresfcn: preallocates a structure with the appropriate fields to
% extract the results of an ode simulation. It is expected that the ode
% function will return a cell array of strings with the appropriate field
% names when called with no arguments.
%
% Syntax
%
% [results, fn, nOut] = prallocresfcn(T, Y, odeevfcn, odeargs)
%
% Input
%
% sol - solution structure as produced by the Matlab ode solvers acting on
%   the solution function.
%
% odeevfcn - function handle or string containing the name of the ode
%   evaluation function used to generate the results. This function must
%   return a cell array of strings containing the names of the variables
%   returned by the function. These will be used to populate a structure
%   containing fields with the same names.
%
% odeargs - (optional) a cell array extra arguments which must be passed to 
%   the ode function.
%
% skip - (optional) a variable to determine if we must extract results for
%   every time in T, or if some can be skipped. The elements
%   1:skip:numel(T) are to be evaluated, and so an appropriate allocation
%   of space will be made for this. If not supplied a value of 1 is used,
%   and space is reserved for every time step. If skip is to be used, but
%   no extra odeargs are required, odeargs can be an empty cell array.
%

    if nargin < 5
        skip = 1;
    end
    
    % get the field names for the output structure by calling the ode
    % function with no inputs. It should be coded to output a cell array of
    % strings containing these if the fourth input argument is a flag with
    % value 0
    fn = feval(odeevfcn, [], [], [], 0);
    
    nOut = numel(fn);
    
    % test the size of the arrays needed for preallocation by calling the
    % ode function
    [testresult{1:nOut}] = feval(odeevfcn, sol.x(1), sol.y(:,1), mc, 1, odeargs{:});
     
    if skip > numel(sol.x) - 1
        rows = 1;
    else
        rows = ceil(numel(sol.x) / skip); 
    end
    
    results.(fn{1}) = zeros(rows, numel(testresult{1}));
    
    % Now preallocate arrays of the correct sizes  
    for i = 2:numel(fn)
        
        results.(fn{i}) = zeros(rows, numel(testresult{i}));
        
    end

end

