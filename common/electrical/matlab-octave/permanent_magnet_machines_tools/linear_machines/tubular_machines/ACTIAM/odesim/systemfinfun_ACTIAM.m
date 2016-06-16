function [design, simoptions] = systemfinfun_ACTIAM (design, simoptions)

    [design, simoptions] = systemfinfun_linear (design, simoptions, @finfun_ACTIAM); 

end