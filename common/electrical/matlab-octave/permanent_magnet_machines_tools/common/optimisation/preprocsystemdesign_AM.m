function design = preprocsystemdesign_AM(design, simoptions, activecoilsperphase)
% processes some common aspects of electrical machine designs in
% preparation for evaluation by a GA

    % set the grid inductance to phase inductance ratio to 0 if not
    % supplied
    design = setfieldifabsent(design, 'LgVLc', 0);
    
    % make minimum possible wire diameter 0.5 mm
    simoptions = setfieldifabsent(simoptions, 'MinWireDiameter',  0.5/1000);
    
    design.Dc = sqrt(4 * (design.Hc * design.Wc * design.CoilFillFactor * design.DcAreaFac) / pi);
    
    if design.Dc < simoptions.MinWireDiameter
        design.Dc = simoptions.MinWireDiameter;
    end
    
    % determine an appropriate coil configuration
    [design.Branches, design.CoilsPerBranch] = branchfac2coilconfig_AM(activecoilsperphase, design.BranchFac);

end

