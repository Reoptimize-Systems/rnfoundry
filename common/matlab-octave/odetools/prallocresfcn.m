function [results, fn, nOut] = prallocresfcn(T, Y, odeevfcn, odeargs, skip)
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
% T - vector of time points as produced by the ode solver acting on the
%   solution function.
%
% Y - matlrix of values of the solution variables at each time step in T.
%   Each column of Y should be a solution variable.
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
    % strings containing these if no input arguments are supplied
    fn = feval(odeevfcn);
    
    nOut = numel(fn); 
    
    % test the size of the arrays needed for preallocation by calling the
    % ode function
    [testresult{1:nOut}] = feval(odeevfcn, T(1), Y(1,:)', odeargs{:});
     
    if skip > numel(T) - 1
        rows = 1;
    else
        rows = ceil(numel(T) / skip); 
    end
    
    results.(fn{1}) = zeros(rows, numel(testresult{1}));
    
    % Now preallocate arrays of the correct sizes  
    for i = 2:numel(fn)
        
        results.(fn{i}) = zeros(rows, numel(testresult{i}));
        
    end

end
