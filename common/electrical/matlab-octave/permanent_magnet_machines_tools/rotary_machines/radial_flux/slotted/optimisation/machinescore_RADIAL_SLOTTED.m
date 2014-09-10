function [score, design, simoptions] = machinescore_RADIAL_SLOTTED(design, simoptions)
 
    [score, design, simoptions] = machinescore_ROTARY(design, simoptions);
    
    % add penalty for cogging torque
    [score, design, simoptions] = ...
        addpenalty_AM(design, simoptions, 'upper', 'CoggingTorquePeak', design.OptimInfo.BaseScore, score);
    
end