function fillfactor = fillfactorcalc(area, coilturns, dc)

    fulldc = conductord2wired(dc);
    
    Ac = pi * (fulldc/2)^2;
    
    fillfactor = area / (coilturns * Ac);
    
end