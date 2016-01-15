function results = oderesults(T, Y, odeevfcn, odeargs, skip, skipfields)
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
% T, Y - these are the output produced by the ODE solver using the solution
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
    [results, fn, nOut] = prallocresfcn(T, Y, odeevfcn, odeargs, skip);
    
    if skip ~= 1
        % store the times the results values occur at if we are skipping
        % some in the results structure
        results.Tskip = T(1:skip:length(T));
    end
    
    k = 0;
    
    outArgs = cell(1, nOut);
    
    % now recalculate the values and put them in the results structure
    for i = 1:skip:length(T)

        k = k + 1;
        
        [outArgs{1:nOut}] = feval(odeevfcn, T(i), Y(i,:)', odeargs{:});
        
        for j = 1:nOut
            if ~ismember (fn{j}, skipfields)
                results.(fn{j})(k,:) = outArgs{j};
            end
        end
    
    end
    
end