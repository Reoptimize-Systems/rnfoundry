function [score, design, simoptions] = machinescore_ROTARY(design, simoptions)
% determines a score for a rotary type electrical machine
%
% Syntax
%
% [score, design, simoptions] = machinescore_ROTARY(design, simoptions)
%


    design = costestimate_AM(design, simoptions);
     
    % evalaute the basic score for modifiaction by the penalty functions,
    % it is expected the basic scoring function will be in the field
    % basescorefcn in the simoptions structure. By default, the function
    % finfun_AM puts the function 'costscore_AM.m' in this field if it is
    % not supplied
    [score, design] = feval(simoptions.basescorefcn, design, simoptions);
    
    design.OptimInfo.BaseScore = score;
    
    [score, design, simoptions]= electricalpenalties_AM(design, simoptions, score);

    [score, design] = masspenalties_ROTARY(design, simoptions, score);

    [score, design, simoptions] = structuralpenalties_AM(design, simoptions, score);
    
    [score, design, simoptions] = systempenalties_ROTARY(design, simoptions, score);
    
end



