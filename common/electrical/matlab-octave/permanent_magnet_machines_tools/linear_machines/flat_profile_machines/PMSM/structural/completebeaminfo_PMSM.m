function BeamInfo = completebeaminfo_PMSM(design, options)

    design.NInnerPoles = design.Poles(1);
    design.NOuterPoles = design.Poles(2);
    
    BeamInfo.StructSides = design.sides;
    
    % get the overall dimensions of the translator frame
    BeamInfo.width = design.ls * options.alphab;
    BeamInfo.depth = 2 * (design.hbf + design.hm + design.g + design.ht + design.hba);
    BeamInfo.height = (design.NOuterPoles .* design.PoleWidth) - design.PoleWidth;
    
    BeamInfo = completebeaminfo_FM(design, options, BeamInfo);
    
end