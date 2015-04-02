function cellstrs = design2pydict_RADIAL (design)

    cellstrs = design2pydict_ROTARY (design);
    
    % create a function for copy-paste convenience
    % create a function for copy-paste convenience
    spm = @(var, val) sprintf ('structdims["Field"]["%s"] = %10.8g * 1000', var, val);
    spa = @(var, val) sprintf ('structdims["Field"]["%s"] = %10.8g', var, val);
    
    cellstrs = [cellstrs;
pytools.trynesteddict('structdims', {'Field', 'Rbi'}, sprintf('%10.8g * 1000', design.Rbi)); 
{ ...
 spm('Rbo', design.Rbo);
[spm('Rmi', design.Rmi), ' * 1.0'];
 spm('Rmo', design.Rmo);
}];


end