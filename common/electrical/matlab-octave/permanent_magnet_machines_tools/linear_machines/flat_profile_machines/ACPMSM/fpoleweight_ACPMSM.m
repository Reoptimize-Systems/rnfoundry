function weight = fpoleweight_ACPMSM(design, evaloptions)
% fpoleweight_PMSM: caculates the weight of a single pole of the field of a
% linear air-cored permanent magnet synchronous machine

    % back iron mass
    weight = design.Taup * design.dbi * design.ls  * evaloptions.FieldIronDensity;
    
    % magnet mass
    weight = weight + (design.lm * design.ls * design.bp * evaloptions.MagnetDensity);
    
    weight = weight * 9.81;

end