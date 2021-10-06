function CoilFillFactor = fillfactorcalc(area, coilturns, dc)

    fulldc = conductord2wired(dc);
    
    Ac = pi * (fulldc/2)^2;
    
    CoilFillFactor = area / (coilturns * Ac);
    
end