function THD = emfthd_AM(slm_fluxlinkage)
% calculates the total harmonic distotion of the emf waveform produced
% by an electrical machine
%
% Syntax
%
% THD = emfthd_AM (slm_fluxlinkage)
%
% Inputs
%
%   slm_fluxlinkage - a periodic slm object fitted to the flux linkage
%       waveform
%

%     % get the fundumental frequency and its harmonics from the slm object
%     % fitted to the flux lnkage with displacement
%     [fundamental_power, harmonics_power] = decomposeemf(slm_fluxlinkage);
% 
%     % Calculate % THD as a ratio of square roots of the total harmonic
%     % power to the fundumental power
%     THD = 100 * sqrt(harmonics_power) / sqrt(fundamental_power);

    % select the number of periods over which to sample
    periods = 20;
    % choose the sampling frequency, there will be Fs samples per period
    Fs = 2048;
    % calculate the number of sample points to use
    N = periods * Fs;
    
    % generate the vector of sample positions
    x = linspace(0, slm_fluxlinkage.period * periods, N);

    % get the voltage waveform at the sample points
    normdphidx = periodicslmeval(x, slm_fluxlinkage, 1);
    
    Nharmonics = 10;
    
    if isoctave () || ~exist ('thd', 'file')
        [~, ~, ~, ~, THDdb] = prettyFFT(normdphidx,Fs,Nharmonics,true,true);
    else
        THDdb = thd (normdphidx, Fs, Nharmonics);
    end
    
    THD = 100 * (10^( THDdb/20 ));

end
