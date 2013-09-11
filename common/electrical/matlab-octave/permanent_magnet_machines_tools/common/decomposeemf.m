function [fundamental_power, harmonics_power, f, harmonics, fundamental_bin] = decomposeemf(slm_fluxlinkage)
% decomposes the derivative of the flux linkage waveform with respect to
% displacement into it fundumental and harmonics
%
% Syntax
%
% [fundamental_power, harmonics_power] = decomposeemf(slm_fluxlinkage)
%
% Inputs
%
%   slm_fluxlinkage - a periodic slm object fitted to the flux linkage
%       waveform
% 

    % select the number of periods over which to sample
    periods = 1;
    % choose the sampling frequency, there will be Fs samples per period
    Fs = 2000;
    % calculate the number of sample points to use
    N = periods * Fs;
    % calculate the bin spacing in Hz
    bin_spacing = Fs / N;
    
    % generate the vector of sample positions
    x = linspace(0, slm_fluxlinkage.period * periods, N);

    % get the voltage waveform at the sample points
    normdphidx = periodicslmeval(x, slm_fluxlinkage, 1);

    % Use fft to find harmonics
    harmonics = abs(real(fft(normdphidx))).^2;

    % remove higher harmonics
    harmonics(101:end) = [];    

    % get the fundumental component (first bin is DC component, the
    % fundumental component is always at a frequency of 1 Hz, per the
    % generation of samples above)
    fundamental_bin = round(1/bin_spacing) + 1;
    n_overtones = 10;
    overtone_bins = round((2:(n_overtones-1))/bin_spacing) + 1;

    overtone_bins(overtone_bins > 100) = [];

    % Find the sum of the rest of the harmonics
    harmonics_power = sum(harmonics(overtone_bins));
    fundamental_power = harmonics(fundamental_bin);
    
    if nargout > 2
        f = bin_spacing * (0:numel(harmonics)-1);
    end

end