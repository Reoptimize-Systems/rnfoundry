function [EMFrms, Irms, Pload] = performanceestimate_ROTARY (design, varargin)

    if isfield (design, 'OmegaPeak')
        % use rpm in last simulation
        options.RPM = omega2rpm (design.OmegaPeak);
    else
        % use rpm which gives 50Hz electrical frequency
        options.RPM = omega2rpm (2 * pi * 50 / (design.Poles/2));
    end
    
    if isfield (design, 'LoadResistance')
        options.LoadResistance = design.LoadResistance;
    else
        options.LoadResistance = 10 * design.PhaseResistance;
    end
    
    if isfield (design, 'LoadInductance')
        options.LoadInductance = design.LoadInductance;
    else
        options.LoadInductance = 0;
    end
    
    if isfield (design, 'LoadCapacitance')
        options.LoadCapacitance = design.LoadCapacitance;
    else
        options.LoadCapacitance = 0;
    end
    
    
    options = parse_pv_pairs (options, varargin);
    

    EMFrms = peakemfest_ROTARY ( design.FluxLinkagePhasePeak, ...
                                 rpm2omega (options.RPM), ...
                                 design.Poles/2 ) / sqrt(2);
                          
                          
	Z = impedance ( design.PhaseResistance(1) + options.LoadResistance, ...
                    design.PhaseInductance(1) + options.LoadInductance, ...
                    options.LoadCapacitance, ...
                    (design.Poles/2) * rpm2omega (options.RPM) / (2*pi) );
                
	Irms = real(EMFrms ./ Z);
    
    Pload = design.Phases .* Irms .* options.LoadResistance;


end