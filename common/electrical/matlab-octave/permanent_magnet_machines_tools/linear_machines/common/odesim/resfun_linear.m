function [results, design] = resfun_linear(T, Y, design, simoptions)
% resfun_linear calculates the results from a linear generator simulation
% performed using the standard linear generator simulation functions

    % we now call the generic resfun_AM function to obtain results
    [results, design] = resfun_AM(T, Y, design, simoptions);
    
    % calculate the input power and efficiency if the required information
    % is present
    if all(isfield(results, {'Fpto', 'vT'}))
        
        if ~isfield(results, 'FaddE')
            results.FaddE = 0;
        end
        
        % instantaneous power is force / velocity
        results.Pinput = -(results.Fpto + results.FaddE) .* results.vT;
        
        design.PowerInputMean = contmean(T, results.Pinput);
        
        design.EnergyInputTotal = trapz(T, results.Pinput);
        
        design.Efficiency = design.EnergyLoadTotal / design.EnergyInputTotal;
        
    end

end
