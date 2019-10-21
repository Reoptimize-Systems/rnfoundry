function a = structfield2array (S, fieldname)
% converts one field of an array of structures to an array of that class

    assert (isfield (S, fieldname), '%s is not a field in the input structure', fieldname);
    
    a = arrayfun (@(x) x.(fieldname), S);

end