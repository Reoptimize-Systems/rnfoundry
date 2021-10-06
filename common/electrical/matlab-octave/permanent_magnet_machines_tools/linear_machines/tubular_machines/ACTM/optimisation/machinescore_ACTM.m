function [score, design] = machinescore_ACTM(design, simoptions)
% score an ACTM machine design based on the cost, power output, and other
% penalties
%

    [score, design] = machinescore_TM(design, simoptions);
    
end