function cellstrs = design2pydict_AM (design)

    % create a function for copy-paste convenience
    spm = @(var, val) sprintf ('structdims["%s"] = %.17e * 1000', var, val);
    spa = @(var, val) sprintf ('structdims["%s"] = %.17e', var, val);
    spd = @(var, val) sprintf ('structdims["%s"] = %d', var, val);
    spAm = @(var, val) sprintf ('structdims["Armature"]["%s"] = %.17e * 1000', var, val);
    spAa = @(var, val) sprintf ('structdims["Armature"]["%s"] = %.17e', var, val);
    spFm = @(var, val) sprintf ('structdims["Field"]["%s"] = %.17e * 1000', var, val);
    spFa = @(var, val) sprintf ('structdims["Field"]["%s"] = %.17e', var, val);
    spFd = @(var, val) sprintf ('structdims["Field"]["%s"] = %d', var, val);

cellstrs = { ...
   'structdims = {}';
   'structdims["Armature"] = {}';
   'structdims["Field"] = {}';
spm('ls', design.ls);
spd('Poles', design.Poles);
spd('Qc', design.Qc);
spFa('MagnetSkew', design.MagnetSkew);
spFd('NSkewMagnetsPerPole', design.NSkewMagnetsPerPole);
};


end