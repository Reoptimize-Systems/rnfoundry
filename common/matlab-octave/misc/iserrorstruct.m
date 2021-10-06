function result = iserrorstruct(candidate)

    result = false;

    if isstruct(candidate)
        if isfield(candidate, {'message', 'identifier', 'stack'})
            result = true;
        end
    end

end