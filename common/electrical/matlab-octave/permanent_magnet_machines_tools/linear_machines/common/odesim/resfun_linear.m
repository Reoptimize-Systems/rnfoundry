function [results, design] = resfun_linear(T, Y, design, simoptions)
% resfun_linear calculates the results from a linear generator simulation
% performed using the standard linear generator simulation functions

    % we now call the generic resfun_AM function to obtain results
    [results, design] = resfun_AM(T, Y, design, simoptions);
    
    % calculate the input power and efficiency if the required information
    % is present
    if all(isfield(results, {'Fpto', 'vT'})) || all(isfield(results, {'Feff', 'vE'}))
        
        if ~isfield(results, 'FaddE')
            results.FaddE = 0;
        end
        
        if isfield (results, 'Fpto')
            ptoforcefieldname = 'Fpto';
            velname = 'vT';
        elseif isfield (results, 'Feff')
            ptoforcefieldname = 'Feff';
            velname = 'vE';
        end
        
        % instantaneous power is force / velocity
        results.Pinput = -(results.(ptoforcefieldname) + results.FaddE) .* results.(velname);
        
        design.PowerInputMean = contmean(T, results.Pinput);
        
        design.EnergyInputTotal = trapz(T, results.Pinput);
        
        design.Efficiency = design.EnergyLoadTotal / design.EnergyInputTotal;
        
        [design.FrequencyPeak, design.VelocityPeak] = freqest_linear(design, results);
        
    end

end


function [freqpeak, velpeak] = freqest_linear(design, results)
% estimates the max electrical frequency from the machine velocity

    if nargin > 1
        % get the max rpm of the generator
        if isfield(results, 'vT')

            velpeak = max(abs(results.vT));

        elseif isfield(results, 'vE')

            velpeak = max(abs(results.vE));
        end
    elseif isfield (design, VelocityPeak)

        velpeak = design.VelocityPeak;

    end
    
    % estimate the electrical frequency at the peak velocity
    freqpeak = velpeak / (2*design.PoleWidth);

end