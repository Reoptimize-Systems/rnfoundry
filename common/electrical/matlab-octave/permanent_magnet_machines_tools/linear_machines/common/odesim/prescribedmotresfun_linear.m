function [results, design] = prescribedmotresfun_linear(T, Y, design, simoptions)

    if ~isfield(simoptions, 'ODEPhaseCurrentCol')
        simoptions.ODEPhaseCurrentCol = 1;
    end

    [results, design] = resfun_linear(T, Y, design, simoptions);

    design.xAmax = 0;
    
    design.vRmax = max(abs(results.vT));
    
    if isfield(results, 'FaddEBD')
        results.FLoss = results.FaddEBD(:,3);
        results.PLoss =  abs(results.FLoss .* results.vT);
        design.PowerLossMean = contmean(T, results.PLoss);
    end

end
