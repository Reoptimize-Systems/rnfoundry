function [results, design] = prescribedmotresfun_ROTARY(T, Y, design, simoptions)

    if ~isfield(simoptions, 'ODEPhaseCurrentCol')
        simoptions.ODEPhaseCurrentCol = 1;
    end
    
    [results, design] = resfun_ROTARY(T, Y, design, simoptions);

end