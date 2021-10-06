function [design, simoptions] = prescribedmotfinfun_TORUS_CORELESS(design, simoptions)
    
    [design, simoptions] = prescribedmotfinfun_ROTARY(design, simoptions, @finfun_TORUS_CORELESS);
    
end