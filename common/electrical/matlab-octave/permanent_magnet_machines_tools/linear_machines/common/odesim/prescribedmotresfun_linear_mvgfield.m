function [results, design] = prescribedmotresfun_linear_mvgfield(T, Y, design, simoptions)

    if ~isfield(simoptions, 'ODEPhaseCurrentCol')
        simoptions.ODEPhaseCurrentCol = 1;
    end
    
    [results, design] = resfun_linear(T, Y, design, simoptions);
    
    results.Fs = results.FaddF(:,1);
    
    results.Fsnap = results.FaddF(:,end);
    
%     if numel(results.FaddF) > 2
%         results.FfA = results.FaddF(:,2);
%     end
    
    if numel(results.FaddF) > 2
        results.FLinearDrag = results.FaddF(:,2);
    end
    
    design.xFmax = max(abs(Y(:,design.phases+1)));

%     design = odeelectricalresults(T, Y(:,1), results.EMF(:,1), design, simoptions);

end
