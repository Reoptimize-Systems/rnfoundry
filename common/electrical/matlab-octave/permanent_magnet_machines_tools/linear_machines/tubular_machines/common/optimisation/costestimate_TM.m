function [cost, design] = costestimate_TM(design, simoptions, buoymass)

    if nargin < 3
        buoymass = 0;
    end
    
    design.BuoyCost = buoymass * simoptions.Evaluation.BuoyMassCost;
    
    [design.PowerConverterCost, design.PowerConverterRating] = ...
                     powerconvertercostest(design.EMFPhasePeak, ...
                                           design.IPhasePeak, ...
                                           design.PowerLoadMean, ...
                                           design.L(1,1), ...
                                           design.vRmax / (2 * design.PoleWidth) );

    design.CopperCost = design.CopperMass * simoptions.Evaluation.CopperCost;

    design.MagnetCost = design.MagnetMass * simoptions.Evaluation.MagnetCost;

    design.ArmatureIronCost = design.ArmatureIronMass * simoptions.Evaluation.ArmatureIronCost;
    
    design.FieldIronCost = design.FieldIronMass * simoptions.Evaluation.FieldIronCost;
                    
    design.ShaftCost = design.ShaftMass * simoptions.Evaluation.StructMaterialCost;
    
    % calculate the estimated cost of the machine, not including converter
    design.MachineCost = design.MagnetCost + design.CopperCost + ...
                         design.ArmatureIronCost + design.FieldIronCost + ...
                         design.ShaftCost;
    
    % Finally calculate the combined cost
    cost = (design.MachineCost + design.PowerConverterCost) * simoptions.Evaluation.nmachines + design.BuoyCost;

    design.CostEstimate = cost;
    
end