function [design, simoptions] = prescribedmotfinfun_ACPMSM(design, simoptions)
    
    [design, simoptions] = prescribedmotfinfun_linear(design, simoptions, @finfun_ACPMSM);

end