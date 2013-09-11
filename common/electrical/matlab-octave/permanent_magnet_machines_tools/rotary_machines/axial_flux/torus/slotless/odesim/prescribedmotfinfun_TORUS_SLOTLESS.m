function [design, simoptions] = prescribedmotfinfun_TORUS_SLOTLESS(design, simoptions)
    
    [design, simoptions] = prescribedmotfinfun_ROTARY(design, simoptions, @finfun_TORUS_SLOTLESS);
    
end