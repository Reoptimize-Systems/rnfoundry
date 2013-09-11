function [design, simoptions] = finfun_ROTARY(design, simoptions)

    % do things specific to ROTARY
    
    % assign the pole width to be the angular width of a pole
    design = setfieldifabsent(design, 'PoleWidth', design.thetap);
    
    % set the power poles, used to calculate force and power etc
    design = setfieldifabsent(design, 'PowerPoles', design.NCoilsPerPhase);
    
    % set the number of machine stages to be one if not supplied
    design = setfieldifabsent(design, 'NStages', 1);
    
    % set the core density to be the same as the back iron material density
    % if available, or to 7500 otherwise
    design.CoreMaterialDensity = 7650;
    
    % set the machine initial conditions to zeros if not present
    simoptions = setfieldifabsent(simoptions, 'IC', zeros(1, design.phases));
    
    % complete the circuit properties
    [design, simoptions] = circuitprops_AM(design, simoptions);
    
    % call finfun_AM to finish off the design
    [design, simoptions] = finfun_AM(design, simoptions);

end 