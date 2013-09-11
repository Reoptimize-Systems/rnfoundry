function [results, design] = systemresfun_linear(T, Y, design, simoptions)
% postprocesses results of a linear machine based WEC ode system simulation 
%
% [results, design] = systemresfun_linear(T, Y, design, simoptions)
%
% 


    % Now obtain internally calculated values by recalling the function
    % with the results, first preallocating the arrays. This is necessary
    % as the ode solver used may take steps while choosing step sizes which
    % do not form part of the solution      
    
    % We may not want to recalculate and plot every single solution step,
    % so we use the skip value to allow us to skip every x solution points,
    % e.g. skip = 2 would calculate only every other solution point.
    if ~isfield(simoptions, 'skip')
        simoptions.skip = 1;
    end
    
    if ~isfield(simoptions, 'ODEPhaseCurrentCol')
        simoptions.ODEPhaseCurrentCol = 5;
    end

    [results, design] = resfun_linear(T, Y, design, simoptions);
    
    if isfield(results, 'xE')
        peakxT = max(abs(results.xE));
    elseif isfield(results, 'xT')
        peakxT = max(abs(results.xT));
    end
    
    if isfield(results, 'vE')
        design.vRmax = max(abs(results.vE));
    elseif isfield(results, 'vT')
        design.vRmax = max(abs(results.vE));
    end
    
    % Actual stator length that would be required to achieve the
    % output (assuming stator is longer part)
    design.minLongMemberLength = 2 * peakxT + (design.PowerPoles * design.PoleWidth);

    design.minLongMemberPoles = ceil(design.minLongMemberLength ./ design.PoleWidth);

    design.minLongMemberLength = design.minLongMemberPoles * design.PoleWidth;
    
    % forces due to machine losses
    results.FLoss = results.FaddEBD(:,3);
    
    % the end stop force
    results.Fes = results.FaddEBD(:,4);
    
    % the buoy system frictional forces
    results.FBuoyFric = results.FaddEBD(:,5);

end
