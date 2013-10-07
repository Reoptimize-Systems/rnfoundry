function [weight] = apoleweight_PMSM(design, evaloptions)

    % get the armature iron mass
    armatureironmass = design.sides * ...
        ((3 * design.ht * design.Wt * design.ls) + (design.Wp * design.hba * design.ls)) * evaloptions.ArmatureIronDensity;
    
    wlength = design.CoilTurns * rectcoilmtl(design.ls, design.Wp - design.Ws, design.Ws);
    
    % Determine the copper cost based on the length of wire in a coil
    % multiplied by the number of Phases and armature Poles and the wire
    % cross-sectional area etc.
    coppermass = design.Phases * wlength * (pi* (design.Dc/2)^2) * design.sides * evaloptions.CopperDensity;
    
    weight = (armatureironmass + coppermass) * 9.81;
    
end