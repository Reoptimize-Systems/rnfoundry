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

    [fundamental_power, harmonics_power] = decomposeemf(slm_fluxlinkage);

    % Calculate THD as a ratio of square roots of powers
    THD = 100 * sqrt(harmonics_power) / sqrt(fundamental_power);

end
