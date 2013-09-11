function design = materialmasses_TORUS_ROTOR(design, simoptions)

    % area of a section of a circle is given by 
    % A = r^2 * theta/2, we subtract the inner section from the outer
    % section
    theta = 2 * pi * design.taumm / design.taupm;
    
    design.MagnetVol = design.tm * 0.5 * (design.Rmo^2 * theta - design.Rmi^2 * theta) * design.NStages * 2;

    design.MagnetMass = design.MagnetVol * simoptions.evaloptions.MagnetDensity;
    
    % calculate the structural material volume abd mass if not already
    % present
    if ~any(isfield(design, {'StructMaterialVolume', 'StructMaterialMass'}))
        
        design.StructMaterialVolume = pi * (design.Rbo^2 - design.Rbi^2) ...
                                      * (2*design.tbio + max(0, design.NStages-1)*design.tbii) ...
                                      + pi * (design.Rbi^2 - design.Rs^2) * (design.NStages*design.outermagsep + 2*design.tbio + max(0, design.NStages-1)*design.tbii + 2*design.tsuppb) ...
                                      + design.tsuppb * design.tausupp * (design.Rbo - design.Rbi) * design.NModuleSupports * design.NModules * design.NStages;
        
        design.StructMaterialMass = design.StructMaterialVolume * simoptions.evaloptions.StructMaterialDensity;
        
    end
    
end