function [score, design] = machinescore_ACPMSM(design, simoptions)

    if isfield(simoptions, 'BuoyParameters')
        [cost, design] = costestimate_ACPMSM(design, simoptions, simoptions.BuoySim.BuoyParameters.mass_external);
    else
        [cost, design] = costestimate_ACPMSM(design, simoptions);
    end
    
        % Estimate the mass of the translator
    if isfield(simoptions, 'StatorPoles')
        % Poles(1) is the armature Poles(2) is the field for the PMSM
        if simoptions.StatorPoles == 1
            design.massT = design.ArmatureIronMass + design.CopperMass;
        elseif simoptions.StatorPoles == 2
            design.massT = design.BackIronMass + design.MagnetMass + design.StructuralMass;
        else
            error('simoptions.StatorPoles must be 1 or 2')
        end
    else
        design.massT = design.BackIronMass + design.MagnetMass + design.StructuralMass;
    end
    
    [score, design] = ddsysscorecommon_linear(design, simoptions);
    
    % exceeding max allowed deflection
    if isfield(design, 'StructPenalty')
        
        [design, simoptions] = structuralpenalty_FM(design, simoptions);
        
        score = score + design.maxAllowedDeflectionpenalty;
        
    end
    
    
    
end