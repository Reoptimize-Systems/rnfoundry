function [cost, design] = costestimate_ACPMSM(design, simoptions, buoymass)
% estimates the cost of the linear air-cored permanent magnet synchronous
% machine

    if nargin < 3
        buoymass = 0;
    end
    
    design.BuoyCost = buoymass * simoptions.evaloptions.BuoyMassCost;
    
    [design.PowerConverterCost, design.PowerConverterRating] = ...
                         powerconvertercostest(design.EMFPhasePeak, ...
                                               design.IPhasePeak, ...
                                               design.PowerLoadMean, ...
                                               design.L(1,1), ...
                                               design.vRmax / (2 * design.PoleWidth) );
                                       
    % get the total length of wire in a coil, based on mean turn length and
    % the number of turns
    design.wlength = design.CoilTurns * rectcoilmtl(design.ls, design.bp, (design.Taup - design.bp) / 2);
    
    % Determine the copper cost based on the length of wire in a coil
    % multiplied by the number of phases and armature poles and the wire
    % cross-sectional area etc.
    design.CopperMass = design.phases * design.wlength * (pi* (design.Dc/2)^2) * design.poles(2) * simoptions.evaloptions.CopperDensity;
    design.CopperCost = design.CopperMass * simoptions.evaloptions.CopperCost;
    
    % Calculate the total cost of the magnets
    design.MagnetMass = design.poles(1) * 2 * design.lm * design.ls * design.bp * simoptions.evaloptions.MagnetDensity;
    design.MagnetCost = design.MagnetMass * simoptions.evaloptions.MagnetCost;
    
    % Calculate the cost of the back iron on the field
    design.FieldIronMass = design.poles(1) * 2 * design.Taup * design.ls * design.dbi * simoptions.evaloptions.FieldIronDensity;;
    design.BackIronMass = design.FieldIronMass; % in case of compatibility issues
    design.FieldIronCost = design.FieldIronMass * simoptions.evaloptions.FieldIronCost;
    
    % Calculate the cost of the structure
    if isfield(design, 'OuterStructureBeamVars')
        design.StructuralMass = structvol_ACPMSM(design, simoptions.evaloptions) * simoptions.evaloptions.StructMaterialDensity;
    else
        design.StructuralMass = 0.5 * design.MagnetMass;
    end
    
    design.structureCost = design.StructuralMass * simoptions.evaloptions.StructMaterialCost;
    
    if isfield(design, 'GuideRailIVars')
        % Calculate the cost of the guide rails
        design.guideRailMass = (design.poles(2) .* design.Taup - design.Taup) * CSArea(design.GuideRailIVars, design.GuideRailIMethod) * 4 * simoptions.evaloptions.StructMaterialDensity;

    else

        design.guideRailMass = 0;

    end
    
    design.guideRailCost = design.guideRailMass * simoptions.evaloptions.StructMaterialCost;
    
    % calculate the estimated cost of the machine, not including converter
    design.machineCost = design.MagnetCost ...
                         + design.CopperCost ...
                         + design.FieldIronCost ...
                         + design.structureCost ...
                         + design.guideRailCost;
                     
    % Finally calculate the combined cost
    cost = (design.machineCost + design.PowerConverterCost) * simoptions.NoOfMachines + design.BuoyCost;
    
    design.CostEstimate = cost;

end