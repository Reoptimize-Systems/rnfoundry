function y1c = translatorY1Centroid(Taup, dbi, bt, hs)

    y1c = ((Taup .* (dbi.^2) / 2) + 3.* hs.* bt .* (dbi + hs/2)) ./ (Taup.*dbi + 3.*hs.*bt);

end