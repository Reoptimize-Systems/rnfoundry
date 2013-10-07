function design = materialmasses_TORUS_SLOTLESS(design, simoptions)

    % get the magnet mass (the rotor structural mass should have been
    % calculated by integrating over the FE structure mesh)
    design = materialmasses_TORUS_ROTOR(design, simoptions);
    
    % calculate the copper mass
    design.CopperVol = design.NStages * design.NCoilsPerPhase * design.Phases * design.MTL * design.CoilTurns * (pi * (design.Dc/2)^2);
    
    design.CopperMass = design.CopperVol * simoptions.evaloptions.CopperDensity;

    % calculate the yoke mass
    design.ArmatureIronVol = pi * (design.Rmo^2 - design.Rmi^2) * design.ty * design.NStages;
    design.ArmatureIronMass = design.ArmatureIronVol * simoptions.evaloptions.ArmatureIronDensity;
    
    % set the epoxy mass to zero
    design.EpoxyVol = 0;
    design.EpoxyMass = design.EpoxyVol * simoptions.evaloptions.EpoxyDensity;
    
    % TODO calculate the aluminium support mass
    design.ArmatureSupportMass = 0;
    
    % calculate the component and module masses
    design.RotorMass = design.StructMaterialMass + design.MagnetMass;
    design.RotorModuleMass = design.RotorMass / design.NModules;
    
    design.StatorMass = design.CopperMass ...
                        + design.EpoxyMass ...
                        + design.ArmatureSupportMass ...
                        + design.ArmatureIronMass;
                    
    design.StatorModuleMass = design.StatorMass / design.NModules;
    
end