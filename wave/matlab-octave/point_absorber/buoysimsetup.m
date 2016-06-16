function buoysimoptions = buoysimsetup(buoy, buoysimoptions)
% buoysimsetup: gets the data for a chosen buoy for use in a simulation
%
% Syntax
%
% buoydat = buoysimsetup(buoy, buoydat)
%
% Input
%
%   buoy - either a scalar, string or empty matrix. If empty, it is assumed
%    that the buoydat structure contains the fields HeaveFile, SurgeFile,
%    ExcitationFile and HydroCoeffsFile, containing the appropriate strings
%    naming the file locations relative to the main buoy library directory.
%
%    If a scalar, the appropriate buoy data is loaded from the directory
%    corresponding to this number in the buoy library directory.
%
%    If a string, this is assumed to be the directory of the chosen buoy,
%    relative to the main buoy library directory.
%
%   buoydat - (optional) a structure which will be returned with the
%    appropriate simulation data appended. If not supplied a new structure
%    will be created with the appropriate data.
%
% Output
%
%   buoydat - A structure containing the following fields
%
%     HeaveFile
%     SurgeFile
%     ExcitationFile
%     HydroCoeffsFile
%     BuoyParameters
%
%

% Copyright Richard Crozier & Helen Bailey 2010-2016

    if nargin == 1
        
        if isstruct(buoy)
            
            buoysimoptions = buoy;
        
            if isfield(buoysimoptions, 'buoynum')
                
                buoy = buoysimoptions.buoynum;
                
            elseif isfield(buoysimoptions, 'buoydir')
                
                buoy = buoysimoptions.buoydir;
                
            else
                error('BUOYSIM:incorrectstruct', 'One input argument was supplied which was a structure without the minimum required fields to specify a simulation.')
            end

        else

            if isscalar(buoy)

                buoysimoptions = struct('buoynum', buoy);

            elseif ischar(buoy)

                buoysimoptions = struct('buoydir', buoy);

            elseif isempty(buoy)
                
                error('BUOYSIM:incorrectinput', 'One argument was supplied which was empty.')
                
            end

        end
    end
    
    if ~isfield(buoysimoptions, 'buoylibdir')
        buoysimoptions.buoylibdir = buoylibdir ();
    end

    if ~isempty(buoy)
        
        if ischar(buoy)
            buoysimoptions = buoydatafromdir(buoysimoptions.buoylibdir, buoy, buoysimoptions);
        elseif isnumeric(buoy)
            if buoy ~= -1
                buoysimoptions = buoynum2buoydata(buoysimoptions.buoylibdir, buoy, buoysimoptions);
            end
        else
            error('Buoy choice is not valid');
        end
        
        if buoy ~= -1
            buoysimoptions.HeaveFile = fullfile(buoysimoptions.buoylibdir, buoysimoptions.HeaveFile);
            buoysimoptions.SurgeFile = fullfile(buoysimoptions.buoylibdir, buoysimoptions.SurgeFile);
            buoysimoptions.ExcitationFile = fullfile(buoysimoptions.buoylibdir, buoysimoptions.ExcitationFile);
            buoysimoptions.HydroCoeffsFile = fullfile(buoysimoptions.buoylibdir, buoysimoptions.HydroCoeffsFile);
        end
        
    end
    
    % set up the sea based on the spec
%     buoydat.SeaParameters = defaultseaparamaters(buoydat.SeaParameters, buoydat.BuoyParameters);

    % check the sea parameters are suitible for the buoy
    if max(buoysimoptions.SeaParameters.sigma / (2*pi)) > buoysimoptions.BuoyParameters.maxfreq
       
        warning('Max sea frequency greater than max buoy frequency, sea will be modified (higher frequencies stripped).')
        
        inds = find((buoysimoptions.SeaParameters.sigma / (2*pi)) <= buoysimoptions.BuoyParameters.maxfreq);
        
        buoysimoptions.SeaParameters.sigma = buoysimoptions.SeaParameters.sigma(inds);
        buoysimoptions.SeaParameters.phase = buoysimoptions.SeaParameters.phase(inds);
        buoysimoptions.SeaParameters.amp = buoysimoptions.SeaParameters.amp(inds);
        buoysimoptions.SeaParameters.L = buoysimoptions.SeaParameters.L(inds);
        buoysimoptions.SeaParameters.wave_number = buoysimoptions.SeaParameters.wave_number(inds);
        
        buoysimoptions.SeaParameters.sigma_range(2) = max(buoysimoptions.SeaParameters.sigma);
        
    end
    
    if min(buoysimoptions.SeaParameters.sigma / (2*pi)) < buoysimoptions.BuoyParameters.maxfreq
       
        warning('Min sea frequency less than min buoy frequency, sea will be modified (lower frequencies stripped).')
        
        inds = find((buoysimoptions.SeaParameters.sigma / (2*pi)) >= buoysimoptions.BuoyParameters.minfreq);
        
        buoysimoptions.SeaParameters.sigma = buoysimoptions.SeaParameters.sigma(inds);
        buoysimoptions.SeaParameters.phase = buoysimoptions.SeaParameters.phase(inds);
        buoysimoptions.SeaParameters.amp = buoysimoptions.SeaParameters.amp(inds);
        buoysimoptions.SeaParameters.L = buoysimoptions.SeaParameters.L(inds);
        buoysimoptions.SeaParameters.wave_number = buoysimoptions.SeaParameters.wave_number(inds);
        
        buoysimoptions.SeaParameters.sigma_range(1) = min(buoysimoptions.SeaParameters.sigma);
        
    end

    % Find the excitation forces
    [buoysimoptions.BuoyParameters.heave_excit_force, buoysimoptions.BuoyParameters.surge_excit_force] = ...
                unitampexcitforces(buoysimoptions.SeaParameters.sigma, buoysimoptions.ExcitationFile);
    
    % load up alphas and betas for this cylinder
    temp = load(buoysimoptions.HeaveFile, 'alpha', 'beta');
    buoysimoptions.BuoyParameters.Halpha = temp.alpha;
    buoysimoptions.BuoyParameters.Hbeta = temp.beta;
    
    temp = load(buoysimoptions.SurgeFile, 'alpha', 'beta');
    
    if isstruct(temp)
        buoysimoptions.BuoyParameters.Salpha = temp.alpha;
        buoysimoptions.BuoyParameters.Sbeta = temp.beta;
    else
        buoysimoptions.BuoyParameters.Salpha = temp(:,1);
        buoysimoptions.BuoyParameters.Sbeta = temp(:,2);
    end

%     % Find the number of alpha and betas used
%     num_coefficients = length(buoydat.BuoyParameters.Halpha);
%     num_coefficients_S = length(buoydat.BuoyParameters.Salpha);

    % Load the variables necessary to calculate the added mass for the
    % buoy, it is expected that this will be an ascii text data file
    hydrocoeffs = load(buoysimoptions.HydroCoeffsFile, '-ascii');

    Hadded_mass = [0 0]; % [sigma, added mass]
    Sadded_mass = [0 0]; % [sigma, added mass]

    for n = 1:length(hydrocoeffs)
        
        if hydrocoeffs(n,2) == 3 && hydrocoeffs(n,3) == 3
            Hadded_mass (end+1, 1) = hydrocoeffs(n,1); % wave number
            Hadded_mass (end, 2) = hydrocoeffs(n,4); % dimensionlised added mass
        end
        
        if hydrocoeffs(n,2) == 1 && hydrocoeffs(n,3) == 1
            Sadded_mass (end+1, 1) = hydrocoeffs(n,1); % wave number
            Sadded_mass (end, 2) = hydrocoeffs(n,4); % dimensionlised added mass
        end
        
    end
    
    Hadded_mass(1,:) = [];% Removes zero first entry
    Sadded_mass(1,:) = [];% Removes zero first entry

    Hadded_mass = [(2 * pi ./ Hadded_mass(:,1)) , Hadded_mass(:,2) .* buoysimoptions.BuoyParameters.rho];
    buoysimoptions.BuoyParameters.HM_infinity = Hadded_mass(end);
    
    Sadded_mass = [(2 * pi ./ Sadded_mass(:,1)) , Sadded_mass(:,2) .* buoysimoptions.BuoyParameters.rho];
    buoysimoptions.BuoyParameters.SM_infinity = Sadded_mass(end);
        
    initial_vel_ext = 0.0;

    buoysimoptions.BuoyParameters.velocity(:,1) = [initial_vel_ext, 0];
    
    % -----------------------------------------------------
    
    % set the number of radiation coefficients to use in the sim if it has
    % not already been specified
    if ~isfield(buoysimoptions, 'NRadiationCoefs')
        buoysimoptions.NRadiationCoefs = numel(buoysimoptions.BuoyParameters.Halpha);
    else
        if buoysimoptions.NRadiationCoefs < 1
            error('You must specify at least one radiation force coeficient.');
        elseif buoysimoptions.NRadiationCoefs < 5
            warning('You have specified the use of less than 5 radiation coefficients, the accuracy may not be reliable');
        end
    end

    % construct the initial conditions of the system
    buoysimoptions.ODESim.SolutionComponents.BuoyPositionHeave.InitialConditions = 0;
    buoysimoptions.ODESim.SolutionComponents.BuoyPositionHeave.AbsTol = buoysimoptions.BuoyParameters.draft / 20;
    
    buoysimoptions.ODESim.SolutionComponents.BuoyVelocityHeave.InitialConditions = 0;
    buoysimoptions.ODESim.SolutionComponents.BuoyVelocityHeave.AbsTol = 0.05;
    
    buoysimoptions.ODESim.SolutionComponents.BuoyPositionSurge.InitialConditions = 0;
    buoysimoptions.ODESim.SolutionComponents.BuoyPositionSurge.AbsTol = buoysimoptions.BuoyParameters.a / 20;
    
    buoysimoptions.ODESim.SolutionComponents.BuoyVelocitySurge.InitialConditions = 0;
    buoysimoptions.ODESim.SolutionComponents.BuoyVelocitySurge.AbsTol = 0.05;
    
    % Handle the  the heave and surge radiation force component setup
    buoysimoptions.ODESim.SolutionComponents.BuoyRadiationHeave.InitialConditions = zeros(1, buoysimoptions.NRadiationCoefs);
    buoysimoptions.ODESim.SolutionComponents.BuoyRadiationSurge.InitialConditions = zeros(1, buoysimoptions.NRadiationCoefs);
    % Set the absolute tolerances on the hydrodynamic forces as the
    % force necessary to accelerate the mass of the buoy at a rate of
    % 0.01 m/s^2
    buoysimoptions.ODESim.SolutionComponents.BuoyRadiationHeave.AbsTol = ...
        repmat(0.01 * buoysimoptions.BuoyParameters.mass_external, 1, buoysimoptions.NRadiationCoefs);
    buoysimoptions.ODESim.SolutionComponents.BuoyRadiationSurge.AbsTol = ...
        repmat(0.01 * buoysimoptions.BuoyParameters.mass_external, 1, buoysimoptions.NRadiationCoefs);
        
    % Calculate a suitible end stop spring constant related to the size of
    % the buoy. 
    %
    % The force from the buoy when submerged by the same amount as its
    % draft is given by:
    %
    %   F = pi * r^2 * rho * g * d
    %
    % where d is the draft. The spring constant is given by
    %
    %   ks = F / x
    %
    % I choose to set the spring constant to be equal to the buoyancy force
    % given previously at some fraction, f, of the buoy draft giving
    %
    %   ks = pi * r^2 * rho * g * d / ( f * d )
    %
    % Which simplifies to 
    %
    %   ks = pi * r^2 * rho * g / f 
    %
    f = 0.05;
    buoysimoptions.EndStopks = (buoysimoptions.BuoyParameters.rho * ...
                         buoysimoptions.BuoyParameters.g * ...
                         pi * buoysimoptions.BuoyParameters.a^2) / f;
                     
    
    % Set heave and surge not to be constrained if not set already
	buoysimoptions.SeaParameters = setfieldifabsent(buoysimoptions.SeaParameters, 'ConstrainSurge', 1);
    buoysimoptions.SeaParameters = setfieldifabsent(buoysimoptions.SeaParameters, 'ConstrainHeave', 1);
    
end