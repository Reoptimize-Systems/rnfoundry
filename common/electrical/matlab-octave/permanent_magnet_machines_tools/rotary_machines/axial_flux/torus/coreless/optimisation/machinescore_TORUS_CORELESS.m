function [score, design, simoptions] = machinescore_TORUS_CORELESS(design, simoptions)

%     design = materialmasses_TORUS_CORELESS(design, simoptions);
    
    [score, design, simoptions] = machinescore_ROTARY(design, simoptions);
    
end


