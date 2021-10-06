function [results, design] = prescribedmotresfun_linear_mvgarm(T, Y, design, simoptions)

    if ~isfield(simoptions, 'ODEPhaseCurrentCol')
        simoptions.ODEPhaseCurrentCol = 1;
    end
    
    [results, design] = resfun_linear(T, Y, design, simoptions);
    
%     results.Fs = results.FaddA(:,1);
    
%     results.Fsnap = results.FaddA(:,end);
    
%     if numel(results.FaddA) > 2
%         results.FfA = results.FaddA(:,2);
%     end
 % from forcefcn_linear_mvgarm: 
 
    % FaddEBD = [Fs, FdragA, FEAFy, FLinearDrag, FfT, FdragT, Fes, Floss]
    % spring force
    results.Fs = results.FaddEBD(:,1);
    % fluid drag forces on armature
    results.FdragA = results.FaddEBD(:,2);
    % magnetic coupling force
    results.Fsnap = results.FaddEBD(:,3);
    % artificial linear drag force 
    results.FLinearDrag = results.FaddEBD(:,4);
    
    
%     if numel(results.FaddA) > 2
%         results.FLinearDrag = results.FaddA(:,2);
%     end
%     
%     if numel(results.FaddA) > 3
%         results.FLoss = -results.FaddA(:,4);
%     end
    
    design.xAmax = max(abs(Y(:,design.Phases+1)));
    design.vRmax = max(abs(results.vT - Y(:,design.Phases+2)));

end
