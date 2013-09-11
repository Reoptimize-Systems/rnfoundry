function [design, simoptions] = simfun_RADIAL(design, simoptions)

    % run the simulation for generic axial flux machines
    [design, simoptions] = simfun_ROTARY(design, simoptions);
    
    % now do stuff specific to RADIAL type machines
    
    % by default we will use a single sided external rotor and internal
    % outward facing stator
    design = setfieldifabsent(design, 'StatorType', 'so');
    
end