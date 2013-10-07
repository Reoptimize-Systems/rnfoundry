function [cost, design] = costestimate_ACTIAM(design, options, buoymass)

    if nargin < 3
        buoymass = 0;
    end
    
    design.buoyCost = buoymass * options.BuoyMassCost;
    
    [design.PowerConverterCost, design.PowerConverterRating] = ...
                     powerconvertercostest(design.EMFPhasePeak, ...
                                           design.IPhasePeak, ...
                                           design.PowerLoadMean, ...
                                           design.L(1,1), ...
                                           design.vRmax / (2 * design.PoleWidth) );
                                       
    design.wlength = design.CoilTurns * pi * 2 * (design.Ri + (design.Hc / 2));
    
    % Determine the copper cost based on the length of wire in a coil
    % multiplied by the number of Phases and armature Poles and the wire
    % cross-sectional area etc.
    design.CopperMass = design.Phases * design.wlength * (pi* (design.Dc/2)^2) ...
                        * design.Poles(2) * options.CopperDensity;
    design.CopperCost = design.CopperMass * options.CopperCost;
    
    % Now calculate the total cost of the magnets
    design.MagnetMass = design.Poles(1) * design.Wm * pi * (design.Rm.^2 - design.Rso^2)*options.MagnetDensity;
    design.MagnetCost = design.MagnetMass * options.MagnetCost;
    
    % Calculate the cost of the laminated iron in the armature
    design.ArmatureIronMass = design.Poles(2) * design.Wp * pi * (design.Ra.^2 - design.Ro.^2) * options.ArmatureIronDensity;
    design.ArmatureIronCost = design.ArmatureIronMass * options.ArmatureIronCost;
    
    % Calculate the cost of the back iron on the field, i.e. the steel
    % discs in the translator
    if design.mode == 1 || design.mode == 3
        cavityvol = steelcavityvol_TM(design);
    else
        cavityvol = 0;
    end
    
    design.FieldIronMass = (design.Ws * pi * (design.Rm^2 - design.Rso^2) - cavityvol) ...
                           * design.Poles(1) * options.FieldIronDensity;
    design.FieldIronCost = design.FieldIronMass * options.FieldIronCost;
    
    % Calculate the cost of the shaft
    design.ShaftMass = (((design.supportLengths(1,2) + design.supportLengths(1,2)) * pi * (design.Rso^2 - design.Rsi^2)) ...
                        + design.Poles(1) * design.Wp * pi * (design.Rso^2 - design.Rsi^2))...
                        * options.StructMaterialDensity;
    design.shaftCost = design.ShaftMass * options.StructMaterialCost;
    
%     design.Pr = powerConverterRating(design.EMFPhasePeak,...
%                                      design.EMFPhasePeak ./ ((design.GridResistance + design.PhaseResistance)),...
%                                      2*pi*max(design.vRmax, 2) / (2 * design.Wp),...
%                                      design.L(1,1));
    
%    Pr = 1.1;
%     
%     % calculate the estimated cost of the converter
%     design.converterCost = 79 * (design.Pr * design.PowerLoadMean/ 1e6) * 1000;
%     

    % calculate the estimated cost of the machine, not including converter
    design.machineCost = design.MagnetCost ...
                         + design.CopperCost ...
                         + design.ArmatureIronCost ...
                         + design.FieldIronCost ...
                         + design.shaftCost;
    
    % Finally calculate the combined cost
    cost = design.machineCost + design.PowerConverterCost + design.buoyCost;

    design.CostEstimate = cost;
       
end