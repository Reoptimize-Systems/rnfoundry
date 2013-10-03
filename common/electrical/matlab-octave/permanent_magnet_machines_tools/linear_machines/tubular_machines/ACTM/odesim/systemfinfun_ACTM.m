function [design, simoptions] = systemfinfun_ACTM(design, simoptions)

	[design, simoptions] = systemfinfun_linear(design, simoptions, @finfun_ACTM);
    
end