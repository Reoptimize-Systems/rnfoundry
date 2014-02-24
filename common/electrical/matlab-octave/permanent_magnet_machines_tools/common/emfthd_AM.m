function THD = emfthd_AM(slm_fluxlinkage)
% calculates the total harmonic distotion of the emf waveform produced
% by an electrical machine
%
% Syntax
%
% THD = thd_AM(slm_fluxlinkage)
%
% Inputs
%
%   slm_fluxlinkage - a periodic slm object fitted to the flux linkage
%       waveform
%

    % get the fundumental frequency and its harmonics from the slm object
    % fitted to the flux lnkage with displacement
    [fundamental_power, harmonics_power] = decomposeemf(slm_fluxlinkage);

    % Calculate % THD as a ratio of square roots of the total harmonic
    % power to the fundumental power
    THD = 100 * sqrt(harmonics_power) / sqrt(fundamental_power);

end
