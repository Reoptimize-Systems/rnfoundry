function [design, simoptions] = systemfinfun_ACPMSM_MC(design, simoptions)

    [design, simoptions] = systemfinfun_linear(design, simoptions, @finfun_ACPMSM_MC);               

end