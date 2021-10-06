function [design, simoptions] = prescribedmotfinfun_RADIAL_SLOTTED (design, simoptions)
    
    [design, simoptions] = prescribedmotfinfun_ROTARY (design, simoptions, @finfun_RADIAL_SLOTTED);
    
end