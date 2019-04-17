function SeaParameters = seasetup(varargin)
% Generates a suitible set of amplitude, phase and angular frequencies for
% use by the hydrodynamic simulation which defines the sea conditions.
%
% Syntax
% 
% SeaParameters = seasetup('Parameter', 'Value')
%
% Input
%
%   The inputs to seasetup are sets of Parameter-Value pairs specifying the
%   properties of the desired sea. The possible parameters are detailed
%   below, grouped into methods of specification:
%
%   Specifying a min and max frequency:
%
%   In this case a set of equally spaced frequencies will be used in the
%   specified range
%
%     'MinFreq' - Minimum frequency in Hz for the frequencies. Defaults to
%        0.1*0.6^0.5 / (2*pi) [which is 12.3281e-003) Hz as recommended
%       by Falco
%   
%     'MaxFreq' - Maximum frequency in Hz for the frequencies, defaults to
%       (0.1*(0.6)^0.5 + 2.24) / (2*pi) [which is 0.3688] Hz as recommended
%       by Falco
%
%    'NoOfFrequencies' - If specifying a minimum and maximum frequency this
%      field determines the number of frequencies generated between the
%      limits. Defualts to 1 if not supplied, and 'FrequencySpacing' is not
%      used instead.
%
%    'FrequencySpacing' - If specifying a minimum and maximum frequency
%      this field determines the spacing between frequencies.
%      NoOfFrequencies is used if not supplied.
%
%   Specifying the angular wave frequencies directly:
%
%   In this case a set of angular wave frequencies can be supplied directly
%
%    'Sigmas' - Vector of one or more angular frequencies of incident waves,
%      if specified, you must also supply Phases and Amplitudes for each
%      frequency.
%
%    'Phases' - Vector of one or more phases, one for each value of Sigmas.
%      If not supplied will be pi/2 for every value of Sigmas supplied. 
%
%    'Amplitudes' - Vector of one or more wave amplitudes, one for each
%      value of Sigmas. If not supplied, will be 0.5m for every value of
%      Sigmas supplied. Note the amplitudes are the maximum wave excursion,
%      i.e. the default of 0.5 results in a wave of height +/- 0.5m.
%
%   To generate a pierson-moskowitz sea state
%
%    'PMPeakFreq' - Scalar value of the peak frequency of the sea state
% 
%    'PMScaleFactor' - Scalar factor by which to scale the wave amplitudes
%      created by the spectrum. Defaults to 1 (no scaling).
%
%  Common Options
%
%    'WaterDepth' - Depth of the water, default is 50 m
%
%    'WaterDensity' - Density of the water, default is 1025 kg/m^3
%
%    'g' - Accceleration due to gravity, default is 9.81 ms^-2
%   
% Output
%
%   SeaParameters, a structure containing the some of the members
%   described in the input and also the wavelengths and wave numbers of the
%   generated frequencies.
%
% Examples
%
% % create a single frequency sea with a frequency of 0.35 Hz. The
% amplitude of the wave will be the default, 0.5m
% SeaParameters = seasetup ('Sigmas', 2 * pi * 0.35);
%  


    % Default freq range taken from Falcao (Phase control through load
    % control of oscillating-body wave energy converters with
    % hydraulic PTO system , Ocean Engineering 2008)
    Inputs.MinFreq = 0.1*0.6^0.5 / (2*pi); % 0.0123 Hz
    Inputs.MaxFreq = (0.1*(0.6)^0.5 + 2.24) / (2*pi); % 0.3688 Hz
    Inputs.NoOfFrequencies = 1;
    Inputs.FrequencySpacing = 0;
    Inputs.ForceFreqRange = false;
    
    % wave angular frequencies
    Inputs.Sigmas = 2 * pi * 1/9;
    % wave phases
    Inputs.Phases = repmat(pi / 2, size(Inputs.Sigmas));
    % wave amplitudes (ignored/replaced for PM Spectrum)
    Inputs.Amplitudes = repmat(0.5, size(Inputs.Sigmas)); % 0.5 m
    Inputs.WaterDepth = 50; % 50 m
    Inputs.WaterDensity = 1025; % 1025 kg/m^3
    Inputs.g = 9.81; % 9.81 ms^-2
    Inputs.PMScaleFactor = 1;

    % settings for a PM Spectrum, if -1 we are not using a PM spectrum
    Inputs.PMPeakFreq = -1;
    
    Inputs.PhaseSeed = [];

    Inputs = parse_pv_pairs(Inputs, varargin);
        
    % restrict the range to that specified by Falco
%     Inputs.MinFreq = max(Inputs.MinFreq, 0.1*0.6^0.5 / (2*pi));
%     Inputs.MaxFreq = min(Inputs.MaxFreq, (0.1*(0.6)^0.5 + 2.24) / (2*pi));

    if Inputs.MaxFreq > (0.1*(0.6)^0.5 + 2.24) / (2*pi)
        warning('SEASETUP:bigmaxfreq', 'Max frequency is greater than that reccomended by Falco');
    end
    
    if Inputs.MinFreq < 0.1*0.6^0.5 / (2*pi)
        warning('SEASETUP:smallminfreq', 'Min frequency is less than that reccomended by Falco');
    end
    
    if Inputs.FrequencySpacing ~= 0
        if Inputs.FrequencySpacing <= (Inputs.MaxFreq - Inputs.MinFreq)
            
            Inputs.NoOfFrequencies = ceil((Inputs.MaxFreq - Inputs.MinFreq) / Inputs.FrequencySpacing);
        
        else
            error('The desired frequency spacing is greater than the specified frequency range.')
        end
    end
    
    if ~samesize(Inputs.Amplitudes, Inputs.Phases)
        
        error(['If supplying phases for the frequencies you must ', ...
               'also supply amplitudes and vice-versa, and the same ', ...
               'number of each.'])
        
    end

    if Inputs.PMPeakFreq == -1

        Inputs.NoOfFrequencies = numel(Inputs.Sigmas);
        
        % if necessary create a random set of phases for the waves if
        % not supplied
        if numel(Inputs.Phases) ~= Inputs.NoOfFrequencies
            Inputs.Phases = 2 * pi * rand(1, Inputs.NoOfFrequencies);
        end

    elseif Inputs.PMPeakFreq > 0

        if Inputs.NoOfFrequencies == 1
            Inputs.NoOfFrequencies = 50;
        end

        SeaParameters.pm_peak_freq = Inputs.PMPeakFreq;
        % [amp, freqs, phase, power]
        [freqs, Inputs.Amplitudes, Inputs.Phases] = pmspectrum(Inputs.PMPeakFreq, ...
                                         'frange', [Inputs.MinFreq, Inputs.MaxFreq], ...
                                         'N', Inputs.NoOfFrequencies, ...
                                         'PhaseSeed', Inputs.PhaseSeed, ...
                                         'ForceFreqRange', Inputs.ForceFreqRange);

        Inputs.Sigmas = 2 * pi * freqs;

        % multiply the wave amplitudes by the scaling factor (normally
        % 1) to make artificially larger waves
        Inputs.Amplitudes = Inputs.Amplitudes * Inputs.PMScaleFactor;                                      

    else
        error('Invalid PM spectrum specifications');
    end

    Inputs.MinSigma = Inputs.MinFreq * 2 * pi;
    Inputs.MaxSigma = Inputs.MaxFreq * 2 * pi;
    Inputs.SigmaRange = [Inputs.MinSigma, Inputs.MaxSigma];

    if any((Inputs.Sigmas > Inputs.MaxSigma) | (Inputs.Sigmas < Inputs.MinSigma))
        error('Desired frequencies are outside the allowed range.')
    end
        
    if numel(Inputs.Sigmas) ~= Inputs.NoOfFrequencies

        if numel(Inputs.Sigmas) ~= 1
            warning('The number of frequencies did not match the number of supplied frequencies')
        end

        % generate appropriate angular frequencies for the sea
        Inputs.Sigmas = linspace(Inputs.MinSigma, Inputs.MaxSigma, Inputs.NoOfFrequencies);

    end
    
    SeaParameters.phase = Inputs.Phases;
    SeaParameters.sigma = Inputs.Sigmas;
    SeaParameters.amp = Inputs.Amplitudes;
    SeaParameters.sigma_range = Inputs.SigmaRange;
    SeaParameters.water_depth = Inputs.WaterDepth;
    
    % calculate the wavelengths of the different frequencies
    SeaParameters.L = seawavelength(SeaParameters.sigma ./ (2*pi), SeaParameters.water_depth, Inputs.g);
    
    % from this determine the wavenumber
    SeaParameters.wave_number = (2*pi) ./ SeaParameters.L;
    
    % Set the water density
    SeaParameters.rho = Inputs.WaterDensity;

end

