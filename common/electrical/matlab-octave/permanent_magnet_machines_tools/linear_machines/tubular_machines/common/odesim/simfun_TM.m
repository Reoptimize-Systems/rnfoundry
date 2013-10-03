function [design, simoptions] = simfun_TM(design, simoptions)

    design.MTL = 2 * pi * mean([design.Ro, design.Ri]); 

    [design, simoptions] = simfun_linear(design, simoptions);

end