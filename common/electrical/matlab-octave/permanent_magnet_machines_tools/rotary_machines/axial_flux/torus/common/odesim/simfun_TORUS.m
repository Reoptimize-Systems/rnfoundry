function [design, simoptions] = simfun_TORUS(design, simoptions)

    % run the simulation for generic axial flux machines
    [design, simoptions] = simfun_AF(design, simoptions);
    
    % now do stuff specific to TORUS type machines (like get forces?)
        
    if design.NStages > 1
        simoptions.ndrawnstages = 2;
    else
        simoptions.ndrawnstages = 1;
    end
    
end