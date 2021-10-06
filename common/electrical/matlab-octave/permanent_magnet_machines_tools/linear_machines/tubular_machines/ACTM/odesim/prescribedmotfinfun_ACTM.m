function [design, simoptions] = prescribedmotfinfun_ACTM(design, simoptions)
    
    [design, simoptions] = prescribedmotfinfun_linear(design, simoptions, @finfun_ACTM);

end