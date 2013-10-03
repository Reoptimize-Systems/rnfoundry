function [cost, design] = costestimate_TM(design, evaloptions, buoymass)

    if nargin < 3
        buoymass = 0;
    end
    
    design.BuoyCost = buoymass * evaloptions.BuoyMassCost;
    
    [design.PowerConverterCost, design.PowerConverterRating] = ...
                     powerconvertercostest(design.EMFPhasePeak, ...
                                           design.IPhasePeak, ...
                                           design.PowerLoadMean, ...
                                           design.L(1,1), ...
                                           design.vRmax / (2 * design.PoleWidth) );
    
    %design.wlength = wirelength_TM(design.RoVRm, design.CoilTurns, design.Rm, design.g);
    design.wlength = design.CoilTurns * pi * 2 * (design.Ri + (design.Hc / 2));
    
    % Determine the copper cost based on the length of wire in a coil
    % multiplied by the number of phases and armature poles and the wire
    % cross-sectional area etc.
    design.CopperMass = design.phases * design.wlength * (pi* (design.Dc/2)^2) * ...
                        design.poles(2) * evaloptions.CopperDensity;
    design.CopperCost = design.CopperMass * evaloptions.CopperCost;
    
    % Now calculate the total cost of the magnets
    design.MagnetMass = design.poles(1) * design.Wm * pi * (design.Rm.^2 - design.Rso^2) * ...
                        evaloptions.MagnetDensity;
    design.MagnetCost = design.MagnetMass * evaloptions.MagnetCost;
    
    % Calculate the cost of the laminated iron in the armature
    design.ArmatureIronMass = design.poles(2) * design.Wp * pi * (design.Ra.^2 - design.Ro.^2) * ...
                              evaloptions.ArmatureIronDensity;
    design.ArmatureIronCost = design.ArmatureIronMass * evaloptions.ArmatureIronCost;
    
    % Calculate the cost of the back iron on the field, i.e. the steel
    % discs in the translator
    if design.mode == 1 || design.mode == 3
        cavityvol = steelcavityvol_TM(design);
    else
        cavityvol = 0;
    end
    
    design.FieldIronMass = (design.Ws * pi * (design.Rm^2 - design.Rso^2) - cavityvol) * ...
                           design.poles(1) * evaloptions.FieldIronDensity;
    design.FieldIronCost = design.FieldIronMass * evaloptions.FieldIronCost;
                    
    design.ShaftCost = design.ShaftMass * evaloptions.StructMaterialCost;
    
    % calculate the estimated cost of the machine, not including converter
    design.MachineCost = design.MagnetCost + design.CopperCost + ...
                         design.ArmatureIronCost + design.FieldIronCost + ...
                         design.ShaftCost;
    
    % Finally calculate the combined cost
    cost = (design.MachineCost + design.PowerConverterCost) * evaloptions.nmachines + design.BuoyCost;

    design.CostEstimate = cost;
    
end