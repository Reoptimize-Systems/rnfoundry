function [design, simoptions] = prescribedmotfinfun_PMSM(design, simoptions)
    
    [design, simoptions] = prescribedmotfinfun_linear(design, simoptions, @finfun_PMSM);

end