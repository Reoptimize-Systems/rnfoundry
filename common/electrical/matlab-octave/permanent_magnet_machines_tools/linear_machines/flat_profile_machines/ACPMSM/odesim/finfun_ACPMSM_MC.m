function [design, simoptions] = finfun_ACPMSM_MC(design, simoptions)

    if ~isfield(simoptions, 'backIronDensity')
        simoptions.backIronDensity = 7500;
    end
    
    if ~isfield(simoptions, 'magnetDensity')
        simoptions.magnetDensity = 7600;
    end
    
    
    [design, simoptions] = finfun_ACPMSM(design, simoptions);
    
    % back iron mass
    ACPMSMFieldPoleMass = design.Taup * design.dbi * design.ls  * simoptions.backIronDensity;
    
    % magnet mass
    ACPMSMFieldPoleMass = ACPMSMFieldPoleMass + (design.lm * design.ls * design.bp * simoptions.magnetDensity);
    
    ACPMSMFieldMass = ACPMSMFieldPoleMass * design.Poles(1) * simoptions.NoOfMachines;
    
    % total mass 
    design.massF = ACPMSMFieldMass + design.MagCouple.FieldMass;
    
    design.xSC = abs(sin(design.AngleFromHorizontal)) * (design.massF * 9.81) / design.ks;
    
%     design.weightF = abs(sin(design.AngleFromHorizontal)) * design.massF * 9.81;
  
end