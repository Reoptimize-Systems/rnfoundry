function weight = apoleweight_ACPMSM(design, evaloptions)
% calculates the weight of a single pole of the field of a linear air-cored
% permanent magnet synchronous machine

    % copper mass
    
    % get the total length of wire in a coil, based on mean turn length and
    % the number of turns
    wlength = design.CoilTurns * rectcoilmtl(design.ls, design.bp, (design.Taup - design.bp) / 2);
    
    % Determine the copper weight based on the length of wire in a coil
    % multiplied by the number of Phases and armature Poles and the wire
    % cross-sectional area etc.
    weight = design.Phases * wlength * (pi* (design.Dc/2)^2) * evaloptions.CopperDensity;
    
    weight = weight * 9.81;

end