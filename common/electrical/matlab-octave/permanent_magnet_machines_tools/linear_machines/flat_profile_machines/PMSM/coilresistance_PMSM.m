function R = coilresistance_PMSM(design)

    tlength = rectcoilmtl(design.ls, design.Wp - design.Wc, design.Wc);
    
    R = design.CoilTurns .* design.WireResistivityBase .* tlength ./ (pi * (design.Dc./2).^2);

end
