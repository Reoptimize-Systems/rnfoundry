function [design, simoptions] = systemfinfun_PMSM(design, simoptions)

    [design, simoptions] = systemfinfun_linear(design, simoptions, @finfun_PMSM); 

end