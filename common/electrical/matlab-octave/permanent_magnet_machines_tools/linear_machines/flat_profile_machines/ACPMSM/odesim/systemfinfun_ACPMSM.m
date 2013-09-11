function [design, simoptions] = systemfinfun_ACPMSM(design, simoptions)

    [design, simoptions] = systemfinfun_linear(design, simoptions, @finfun_ACPMSM);             

end