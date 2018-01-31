function [cost, design] = costestimate_PMSM(design, simoptions, buoymass)
% costestimate_PMSM, estimates the cost of the linear permanent magnet
% synchronous machine

    if nargin < 3
        buoymass = 0;
    end
    
    % Estimate the cost of the power converter required
    [design.PowerConverterCost, design.PowerConverterRating] = ...
                     powerconvertercostest(design.EMFPhasePeak, ...
                                           design.IPhasePeak, ...
                                           design.PowerLoadMean, ...
                                           design.L(1,1), ...
                                           design.vRmax / (2 * design.PoleWidth) );
    
    design.BuoyCost = buoymass * simoptions.Evaluation.BuoyMassCost;
    
    % design.wlength = GeWcireLength(design.RoVRm, design.CoilTurns, design.Rm, design.g);
    
    design.wlength = design.CoilTurns * rectcoilmtl(design.ls, design.Wp - design.Wc, design.Wc);
    
    % Determine the copper cost based on the length of wire in a coil
    % multiplied by the number of Phases and armature Poles and the wire
    % cross-sectional area etc.
    design.CopperMass = design.Phases * design.wlength * (pi* (design.Dc/2)^2) ...
                        * design.Poles(1) * design.sides * simoptions.Evaluation.CopperDensity;
    design.CopperCost = design.CopperMass * simoptions.Evaluation.CopperCost;
    
    % Calculate the total cost of the magnets
    design.MagnetMass = design.Poles(2) * design.sides * design.hm * design.ls ...
                        * design.Wp * simoptions.Evaluation.MagnetDensity;
    design.MagnetCost =  design.MagnetMass * simoptions.Evaluation.MagnetCost;
    
    % Calculate the cost of the laminated iron (teeth and yoke) in the armature
    design.ArmatureIronMass = design.Poles(1) * design.sides * ((3 * design.ht * design.Wt * design.ls) ...
                              + (design.Wp * design.hba * design.ls)) * simoptions.Evaluation.ArmatureIronDensity;
    design.ArmatureIronCost = design.ArmatureIronMass * simoptions.Evaluation.ArmatureIronCost;
    
    % Calculate the cost of the back iron on the field
    design.BackIronMass = design.Poles(2) * design.sides * design.Wp * design.ls ...
                          * design.hbf * simoptions.Evaluation.FieldIronDensity;
    design.BackIronCost = design.BackIronMass * simoptions.Evaluation.FieldIronCost;
    
    if every(isfield(design, {'GuideRailIVars', 'GuideRailIMethod'}))
        
        % Calculate the cost of the support structure
        design.StructuralMass = structvol_PMSM(design, simoptions.Evaluation) * simoptions.Evaluation.StructMaterialDensity;

        design.StructuralCost = design.StructuralMass * simoptions.Evaluation.StructMaterialCost;

        % Calculate the cost of the guide rails
        design.guideRailMass = (design.Poles(1) .* design.Wp - design.Wp) * CSArea(design.GuideRailIVars, design.GuideRailIMethod) * 2 * design.sides * simoptions.Evaluation.StructMaterialDensity;

        design.guideRailCost = design.guideRailMass * simoptions.Evaluation.StructMaterialCost;
        
    else
        
        design.StructuralCost = 0;
        
        design.StructuralMass = 0;
        
        design.guideRailCost = 0;
        
        design.guideRailMass = 0;
        
    end
    
    design.massT = design.StructuralMass + design.BackIronMass + design.MagnetMass;
    
    % calculate the estimated cost of the machine, not including converter
    design.machineCost = design.MagnetCost ...
                         + design.CopperCost ...
                         + design.ArmatureIronCost ...
                         + design.BackIronCost ...
                         + design.StructuralCost ...
                         + design.guideRailCost;
    
    % Finally calculate the combined cost
    cost = (design.machineCost + design.PowerConverterCost) * simoptions.NoOfMachines + design.BuoyCost;

    design.CostEstimate = cost;

end