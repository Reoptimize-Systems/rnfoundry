function BeamInfo = completebeaminfo_ACPMSM(design, options)

    design.NInnerPoles = design.poles(2);
    design.NOuterPoles = design.poles(1);
    
    BeamInfo.StructSides = 2;
    
    % get the overall dimensions of the translator frame
    BeamInfo.width = design.ls * options.alphab;
    BeamInfo.depth = 2 * (design.dg + design.lm + design.dbi);
    BeamInfo.height = (design.NOuterPoles .* design.PoleWidth) - design.PoleWidth;
    
    BeamInfo = completebeaminfo_FM(design, options, BeamInfo);
    
end