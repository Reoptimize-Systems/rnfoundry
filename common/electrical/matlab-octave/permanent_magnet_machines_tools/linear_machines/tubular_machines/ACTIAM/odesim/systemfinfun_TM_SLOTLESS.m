function [design, simoptions] = systemfinfun_TM_SLOTLESS (design, simoptions)

	[design, simoptions] = systemfinfun_linear(design, simoptions, @finfun_TM_SLOTLESS );
    
end