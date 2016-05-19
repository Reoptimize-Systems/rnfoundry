function design = coilresistance_RADIAL_SLOTTED (design)

    if  ~isfield (design, 'MTL')
        extralen = design.thetas * design.Rcm / 2;

        design.MTL = rectcoilmtl ( design.ls, ...
                                   design.yd * design.thetas * design.Rcm + extralen, ...
                                   mean (design.thetac) * design.Rcm );
    end
    
    
    design.CoilResistance = wireresistancedc ('round', design.Dc, design.MTL*design.CoilTurns);
    
end