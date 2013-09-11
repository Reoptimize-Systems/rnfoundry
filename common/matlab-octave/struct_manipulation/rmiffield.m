function struct = rmiffield(struct, fields)
% rmiffield: removed a field or several fields from a structure if they
% exists, or otherwise do nothing with no notification

    if ischar(fields)
        fields = {fields};
    end

    struct = rmfield(struct, fields(isfield(struct, fields)));

end