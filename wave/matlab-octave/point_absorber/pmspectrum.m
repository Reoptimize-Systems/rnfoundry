function [freqs, amps, phase, power] = pmspectrum(fpeak, varargin)
% pmspectrum: generates a set of frequencies, phases and wave amplitudes
% corresponding to a pierson-moskowitz sea state as specified by a given
% peak frequency.
%
% Syntax
%
% [amp, freqs, phase, power] = pmspectrum(fpeak, 'param', 'value')
%
%
% Inputs
%
% fpeak - a scalar value of the peak frequency of the desired spectrum in
%   Hz
%
% The following optional parameter-value pairs may then also be supplied:
%
%  'frange' - (1 x 2) vector containing the minimum and maximum values
%     of the desired frequency range of the spectrum. The default range is
%     0.0123 Hz to 0.369 Hz taken from Falcao (Phase control through load
%     control of oscillating-body wave energy converters with hydraulic PTO
%     system , Ocean Engineering 2008)
%
%  'N' - a scalar value of the number of frequencies in the output
%     spectrum, the default is 100
%
%  'g' - a scalar value of the local acceleration due to gravity, the
%     default value is 9.81 ms^-2
%
% Output
%
% freqs - a set of N frequencies linearly spaced between the maximum and
%         minimum values of the desired frequency range.
%
% amps - A set of N wave amplitudes one for for each frequency in freqs
%  
% phase - A set of N random phases offsets, one for each frequency, these
%         are derived from a uniform random distribution.
% 
% power - the spectrum power for each frequency
%

    % Set up the default values for the optional parameters.
    
    % Default frequency range taken from Falcao (Phase control through load
    % control of oscillating-body wave energy converters with hydraulic PTO
    % system , Ocean Engineering 2008)
    Inputs.frange = [0.1*0.6^0.5, 0.1*(0.6)^0.5 + 2.24] / (2*pi);
    % Earth's gravity by default, change this when using on Mars/Titan etc.
    Inputs.g = 9.81;
    % 100 frequencies by default
    Inputs.N = 100;
    Inputs.ForceFreqRange = false;
    
    Inputs.PhaseSeed = [];
    
    % parse the input parameter-value pairs
    Inputs = parse_pv_pairs(Inputs, varargin);

    % create a set of equally spaced frequencies
    freqs = linspace(Inputs.frange(1), Inputs.frange(2), Inputs.N);
    
    % get the interfrequency spacing
    fspacing = freqs(2) - freqs(1);
    
    % Calculate the spectrum power from the frequencies
    power = 8.1e-3 .* Inputs.g.^2 .* ( 2 .* pi ).^-4 .* freqs.^-5 .* exp(-5/4 .* (fpeak ./ freqs).^4);

    % determine the wave amplitudes
    amps = 0.5 .* sqrt(2 .* power .* 2 .* pi .* fspacing);

    % create random phase differences for each frequency
    if isempty (Inputs.PhaseSeed)
        rng('shuffle');
    else
        rng(Inputs.PhaseSeed);
    end
    
    phase = 2 * pi * rand(size(freqs));

end