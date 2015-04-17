function cellstrs = design2pydict_RADIAL_SLOTTED (design)

    cellstrs = design2pydict_RADIAL (design);
    
    % create a function for copy-paste convenience
    spm = @(var, val) sprintf ('structdims["Armature"]["%s"] = %.17e * 1000', var, val);
    spa = @(var, val) sprintf ('structdims["Armature"]["%s"] = %.17e', var, val);

cellstrs = [cellstrs;
pytools.trynesteddict('structdims', {'Armature', 'Rai'}, sprintf('%.17e * 1000', design.Rai)); 
{ ...
spm('Rtsg', design.Rtsg);
spm('Rci', design.Rci);
spm('Rcb', design.Rcb);
spm('Ryi', design.Ryi);
spm('Ryo', design.Ryo);
spa('thetas', design.thetas);
spa('thetasg', design.thetasg);
spa('thetacg', design.thetacg);
spa('thetacy', design.thetacy);
}];

end