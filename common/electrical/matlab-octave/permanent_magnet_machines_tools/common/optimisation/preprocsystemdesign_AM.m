function design = preprocsystemdesign_AM(design, simoptions, activecoilsperphase)
% processes some common aspects of electrical machine designs in
% preparation for evaluation by a GA

    % set the grid inductance to phase inductance ratio to 0 if not
    % supplied
    design = setfieldifabsent(design, 'LgVLc', 0);
    
    % make minimum possible wire diameter 0.5 mm if not already set
    simoptions = setfieldifabsent(simoptions, 'MinStrandDiameter',  0.5/1000);
    % set infinite max wire diameter if not yet set
    simoptions = setfieldifabsent(simoptions, 'MaxStrandDiameter',  inf);
    
    if simoptions.MaxStrandDiameter < simoptions.MinStrandDiameter
        error ('In simoptions, maximum allowed wire size is less than minimum allowed wire size.')
    end
    
    design.Dc = sqrt(4 * (design.Hc * design.Wc * design.CoilFillFactor * design.DcAreaFac) / pi);
    
    % Calculate what strand diameter gives the exact same copper conductor
    % CS area
    design.WireStrandDiameter = stranddiameter (design.Dc, design.NStrands);
    
    % Check it's not smaller than the minimum possible wire size
    if design.WireStrandDiameter < simoptions.MinStrandDiameter

        design.WireStrandDiameter = simoptions.MinStrandDiameter;

        design.Dc = equivDcfromstranded (design.WireStrandDiameter, design.NStrands);

    end
    
    % check it's not bigger than the biggest allowed wire size
    if design.WireStrandDiameter > simoptions.MaxStrandDiameter

        design.WireStrandDiameter = simoptions.MaxStrandDiameter;

        [design.NStrands, design.WireStrandDiameter] = CoilTurns (circlearea(design.Dc/2), 1.0, design.WireStrandDiameter);

        design.Dc = equivDcfromstranded (design.WireStrandDiameter, design.NStrands);

    end
    
    % determine an appropriate series/parallel coil configuration
    [design.Branches, design.CoilsPerBranch] = branchfac2coilconfig_AM (activecoilsperphase, design.BranchFac);

end

