function S = pvpairs2struct (pv_pairs)
% converts a set of parameter-value pairs to a structure

    npv = length(pv_pairs);
    n = npv/2;

    if n~=floor(n)
      error ('Property/value pairs must come in PAIRS.')
    end
    
    fieldnames = pv_pairs(1:2:end);
    params = pv_pairs(2:2:end);
    
    if ~iscellstr (fieldnames)
        error ('Parameter name was not a string');
    end
    
    for ind = 1:numel (fieldnames)
        S.(fieldnames{ind}) = params{ind};
    end
    
end