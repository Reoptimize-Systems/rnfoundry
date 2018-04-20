function call_resets_recurse (ODESim)
% recursively call the reset fuction for the solution components and any
% nested simulation components if present

    if isfield (ODESim, 'NestedSim')
        nested_sim_recurse (ODESim);
    end
    
    if isfield (ODESim, 'SolutionComponents')
        fnames = fieldnames (ODESim.SolutionComponents);
        for ind = 1:numel(fnames)
            if isfield (ODESim.SolutionComponents.(fnames{ind}), 'ResetFcn')
                feval (ODESim.SolutionComponents.(fnames{ind}).ResetFcn);
            end
        end
    end

end