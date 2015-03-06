function [score, design, simoptions] = machinescore_RADIAL_SLOTTED(design, simoptions)
% generates a score for a slotted radial flux permanent magnet machine
% desgn
%
% Syntax
%
% [score, design, simoptions] = machinescore_RADIAL_SLOTTED(design, simoptions)
%
% 
 
    [score, design, simoptions] = machinescore_ROTARY(design, simoptions);
    
    % add penalty for cogging torque
    [score, design, simoptions] = ...
        addpenalty_AM(design, simoptions, 'upper', 'CoggingTorquePeak', design.OptimInfo.BaseScore, score);
    
end