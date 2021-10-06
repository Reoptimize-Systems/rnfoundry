function [design, simoptions] = prescribedmotfinfun_ACTIAM(design, simoptions)

    [design, simoptions] = prescribedmotfinfun_linear(design, simoptions, @finfun_ACTIAM);

end