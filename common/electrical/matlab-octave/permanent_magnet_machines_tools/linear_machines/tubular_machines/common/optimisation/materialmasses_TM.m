function [design, simoptions] = materialmasses_TM(design, simoptions)
% calculates the masses of common components of tubular linear machines
%
% Syntax
%
% [design, simoptions] = materialmasses_TM(design, simoptions)
%
% 

    design.wlength = design.CoilTurns * pi * 2 * (design.Ri + (design.Hc / 2));
    
    % Determine the copper cost based on the length of wire in a coil
    % multiplied by the number of Phases and armature Poles and the wire
    % cross-sectional area etc.
    design.CopperMass = design.Phases * design.wlength * (pi* (design.Dc/2)^2) ...
                        * design.Poles(2) * simoptions.evaloptions.CopperDensity;
    
    % Now calculate the total cost of the magnets
    design.MagnetMass = design.Poles(1) * design.Wm * pi * (design.Rm.^2 - design.Rso^2) * simoptions.evaloptions.MagnetDensity;
    
    % Calculate the cost of the laminated iron in the armature
    design.ArmatureIronMass = design.Poles(2) * design.Wp * pi * (design.Ra.^2 - design.Ro.^2) * simoptions.evaloptions.ArmatureIronDensity;
    
    % Calculate the cost of the back iron on the field, i.e. the steel
    % discs in the translator
    if design.mode == 1 || design.mode == 3
        cavityvol = steelcavityvol_TM(design);
    else
        cavityvol = 0;
    end
    
    design.FieldIronMass = (design.Ws * pi * (design.Rm^2 - design.Rso^2) - cavityvol) ...
                            * design.Poles(1) * simoptions.evaloptions.FieldIronDensity;
    
end