function design = preprocsystemdesign_ROTARY(design, simoptions, activecoilsperphase)
% preprocsystemdesign_ROTARY: processes some common aspects of linear
% machine designs in preparation for evaluation by a GA
    
    design = preprocsystemdesign_AM(design, simoptions, activecoilsperphase);
    
end

