function phi = cricwavephasediffs(L, radius, Nbuoys)

    buoycircangles = linspace(0, 2*pi, Nbuoys+1);
    
    buoycircangles(end) = [];
    
    xdisp = radius .* cos(buoycircangles);
    
    phi = zeros(Nbuoys, numel(L));
    
    for i = 1:numel(xdisp)
        
        phi(i,:) = 2 .* pi .* xdisp(i) ./ L; 
    
    end

end