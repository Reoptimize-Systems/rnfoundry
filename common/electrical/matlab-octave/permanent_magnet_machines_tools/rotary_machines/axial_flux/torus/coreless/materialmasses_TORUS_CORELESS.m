function design = materialmasses_TORUS_CORELESS(design, simoptions)

    % get the magnet mass (the rotor structural mass should have been
    % calculated by integrating over the FE structure mesh)
    design = materialmasses_TORUS_ROTOR(design, simoptions);
    
    % calculate the copper mass
    design.CopperVol = design.NStages * design.NCoilsPerPhase * design.phases * design.MTL * design.CoilTurns * (pi * (design.Dc/2)^2);
    
    design.CopperMass = design.CopperVol * simoptions.evaloptions.CopperDensity;

    % calculate the epoxy mass
    coilthickness = (design.tauco - design.tauci); %/2;
    epoxyregionvol = design.NStages * design.tc * pi * ((design.Rmo + coilthickness)^2 - (design.Rmi - coilthickness)^2);
    design.EpoxyVol = epoxyregionvol - design.CopperVol;
    if design.EpoxyVol < 0
        design.EpoxyVol = 0;
    end
    design.EpoxyMass = design.EpoxyVol * simoptions.evaloptions.EpoxyDensity;
    
    % TODO calculate the aluminium support mass
    design.ArmatureSupportMass = 0;
    
    % calculate the component and module masses
    design.RotorMass = design.StructMaterialMass + design.MagnetMass;
    design.RotorModuleMass = design.RotorMass / design.NModules;
    
    design.StatorMass = design.CopperMass + design.EpoxyMass + design.ArmatureSupportMass;
    design.StatorModuleMass = design.StatorMass / design.NModules;
    
end