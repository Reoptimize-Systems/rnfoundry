function [results, design, summary_time_inds] = resfun_ROTARY(T, Y, design, simoptions)
% post processes results from an ode simulation of a rotary electrical
% machine
%
% Syntax
%
% [results, design] = resfun_ROTARY(T, Y, design, simoptions)
%
%
% Input
%
% 

% Copyright Richard Crozier 2012

    if ~isfield(simoptions, 'ODEPhaseCurrentCol')
        simoptions.ODEPhaseCurrentCol = 1;
    end

    [results, design, summary_time_inds] = resfun_AM(T, Y, design, simoptions);
    
    % calculate the input power and efficiency if the required information
    % is present
    if isfield(results, 'Tqpto')
        
        design.TorquePtoPeak = max(abs(results.Tqpto(summary_time_inds,:)));
            
        if isfield (design, 'Rmm')
            % get the maximum force on the magnets
            %
            % convert torques to forces at the mean magnet radius
            results.Fpto = results.Tqpto ./ design.Rmm;
            %
            % determine the maximum force on the magnets
            design.MaxFpto = design.TorquePtoPeak ./ design.Rmm;
        end
        
        % divide the maximum amplitude of the normalised torque ripple by the
        % mean to get the ripple factor
        design.TorqueRippleFactor = torqueripple(T, results, summary_time_inds);
        
        if isfield (results, 'omegaT')

            if ~isfield(results, 'TqaddE')
                results.TqaddE = 0;
            end

            % instantaneous power is force / velocity
            results.Pinput = -(results.Tqpto + results.TqaddE) .* results.omegaT;

            design.PowerInputMean = contmean(T(summary_time_inds), results.Pinput(summary_time_inds,:));

            design.EnergyInputTotal = trapz(T(summary_time_inds), results.Pinput(summary_time_inds,:));

            if isfield (design, 'EnergyLoadTotal')
                design.Efficiency = design.EnergyLoadTotal / design.EnergyInputTotal;
            end

            design.TorquePtoMean = contmean(T(summary_time_inds), results.Tqpto(summary_time_inds,:));

        end
        
    end

    if isfield(results, 'TqaddEBD')
        
        % get iron losses forces
        results.TqLiron = results.TqaddEBD(:,3);
        results.PLiron =  abs(results.TqLiron .* results.omegaT);
        design.PowerLossIronMean = contmean(T, results.PLiron);
        
        % get winding eddy current losses
        results.TqLeddy = results.TqaddEBD(:,4);
        results.PLeddy =  abs(results.TqLeddy .* results.omegaT);
        design.PowerLossEddyMean = contmean(T, results.PLeddy);
        
        % calculate the mean total non-copper losses
        design.PowerLossMean = contmean(T, results.PLiron + results.PLeddy);
        
    end
    
    [design.FrequencyPeak, design.OmegaPeak] = freqest_ROTARY(design, results, summary_time_inds);

end

function TR = torqueripple(T, results, summary_time_inds)
% calculates a torque ripple factor 
%

    % calculate a torque ripple factor
    
    % normalise the torque values to the tangential velocity
    normripple = results.Tqpto ./ results.omegaT;
    
    % strip values for zero velocity
    normripple(results.omegaT == 0) = 0;
    
    % remove the first 20% of the time to allow for start up conditions
%     normripple = normripple(T > (max(T) - min(T))*0.2);
%     T = T(T > (max(T) - min(T))*0.2);
    
    % divide the maximum amplitude of the normalised torque ripple by the
    % mean to get the ripple factor
    TR = abs((max(normripple(summary_time_inds,:)) - min(normripple(summary_time_inds,:))) / contmean(T(summary_time_inds), normripple(summary_time_inds,:)));

end


function [freqpeak, omegapeak] = freqest_ROTARY(design, results, summary_time_inds)
% estimates the max electrical frequency from the machine velocity

    % get the max rpm of the generator
    if isfield(results, 'vT')

        omegapeak = max(abs(vel2omega(results.vT(summary_time_inds,:), design.Rmm)));

    elseif isfield(results, 'vE')

        omegapeak = max(abs(vel2omega(results.vE(summary_time_inds,:), design.Rmm)));

    elseif isfield(results, 'omegaT')

        omegapeak = max(abs(results.omegaT(summary_time_inds,:)));

	elseif isfield(results, 'omegaE')
        
        omegapeak = max(abs(results.omegaE(summary_time_inds,:)));
        
    end
    
    % get the frequency of rotation
    freqpeak = omegapeak ./ (2*pi);
    
    % estimate the electrical frequency at the peak velocity
    freqpeak = freqpeak .* (design.Poles/2);

end