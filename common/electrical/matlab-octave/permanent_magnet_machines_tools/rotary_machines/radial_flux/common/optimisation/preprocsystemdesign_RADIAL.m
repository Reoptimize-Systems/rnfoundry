function design = preprocsystemdesign_RADIAL(design, simoptions)
% preprocsystemdesign_RADIAL: processes some common aspects of radial flux
% machine designs in preparation for evaluation by a GA
    
    design = preprocsystemdesign_ROTARY(design, simoptions, design.NCoilsPerPhase);
    
end

