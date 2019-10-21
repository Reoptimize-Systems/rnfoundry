function a = structfield2cell (S, fieldname)
% converts one field of an array of structures to a cell array 

    assert (isfield (S, fieldname), '%s is not a field in the input structure', fieldname);
    
    a = arrayfun (@(x) x.(fieldname), S, 'UniformOutput', false);

end