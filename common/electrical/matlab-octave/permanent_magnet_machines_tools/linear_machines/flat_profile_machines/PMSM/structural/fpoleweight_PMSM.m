function weight = fpoleweight_PMSM(design, evaloptions)
% fpoleweight_PMSM: caculates the weight of a single pole of the field of a
% linear permanent magnet synchronous machine

    % back iron mass
    weight = design.Wp * design.hbf * design.ls  * evaloptions.FieldIronDensity;
    
    % magnet mass
    weight = weight + (design.hm * design.ls * design.Wm * evaloptions.MagnetDensity);
    
    weight = weight * 9.81;

end