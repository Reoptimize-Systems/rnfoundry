function [results, design] = prescribedmotresfun_RADIAL_SLOTTED(T, Y, design, simoptions)
% post-processes the results of a prescribed motion ODE simulation of a
% rotary electrical machine
%

    if ~isfield(simoptions, 'ODEPhaseCurrentCol')
        simoptions.ODEPhaseCurrentCol = 1;
    end
    
    [results, design] = resfun_RADIAL_SLOTTED(T, Y, design, simoptions);

end