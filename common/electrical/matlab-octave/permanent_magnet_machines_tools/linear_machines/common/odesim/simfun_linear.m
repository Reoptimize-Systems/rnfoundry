function [design, simoptions] = simfun_linear(design, simoptions)

    [design, simoptions] = simfun_AM(design, simoptions);
    
    design = setfieldifabsent(design, 'AngleFromHorizontal', pi/2);

end