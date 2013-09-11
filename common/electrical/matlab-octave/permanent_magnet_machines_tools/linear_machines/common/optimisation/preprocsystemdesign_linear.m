function design = preprocsystemdesign_linear(design, simoptions, powerpoles)
% preprocsystemdesign_linear: processes some common aspects of linear
% machine designs in preparation for evaluation by a GA

    % set default machine mounting angle
    design.AngleFromHorizontal = 80 * (pi/180);
    
    % set mass of translator to zero initially
    design.massT = 0;
    
    design = preprocsystemdesign_AM(design, simoptions, powerpoles);
    
end

