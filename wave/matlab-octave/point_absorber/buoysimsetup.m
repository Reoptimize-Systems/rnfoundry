function [buoydat, I_H, I_S] = buoysimsetup(buoy, buoydat)
% buoysimsetup: gets the data for the chosen buoy for use in a simulation
%
% Syntax
%
% buoydat = buoysimsetup(buoy, buoydat)
%
% Input
%
%   buoy - either a scalar, string or empty matrix. If empty, it is assumed
%          that the buoydat structure contains the fields HeaveFile,
%          SurgeFile, ExcitationFile and HydroCoeffsFile, containing the
%          appropriate strings naming the file locations relative to the
%          main buoy library directory.
%
%          If a scalar, the appropriate buoy data is loaded from the
%          directory corresponding to this number in the buoy library
%          directory.
%
%          If a string, this is assumed to be the directory of the chosen
%          buoy, relative to the main buoy library directory.
%
%   buoydat - (optional) a structure which will be returned with the
%          appropriate simulation data appended.
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

% Copyright Richard Crozier & Helen Bailey 2010

    if nargin == 1
        
        if isstruct(buoy)
            
            buoydat = buoy;
        
            if isfield(buoydat, 'buoynum')
                
                buoy = buoydat.buoynum;
                
            elseif isfield(buoydat, 'buoydir')
                
                buoy = buoydat.buoydir;
                
            else
                error('BUOYSIM:incorrectstruct', 'One input argument was supplied which was a structure without the minimum required fields to specify a simulation.')
            end

        else

            if isscalar(buoy)

                buoydat = struct('buoynum', buoy);

            elseif ischar(buoy)

                buoydat = struct('buoydir', buoy);

            elseif isempty(buoy)
                
                error('BUOYSIM:incorrectinput', 'One argument was supplied which was empty.')
                
            end

        end
    end
    
    if ~isfield(buoydat, 'buoylibdir')
        buoydat.buoylibdir = getbuoylibdir;
    end

    if ~isempty(buoy)
        
        if ischar(buoy)
            buoydat = buoydatafromdir(buoydat.buoylibdir, buoy, buoydat);
        elseif isnumeric(buoy)
            if buoy ~= -1
                buoydat = buoynum2buoydata(buoydat.buoylibdir, buoy, buoydat);
            end
        else
            error('Buoy choice is not valid');
        end
        
        if buoy ~= -1
            buoydat.HeaveFile = fullfile(buoydat.buoylibdir, buoydat.HeaveFile);
            buoydat.SurgeFile = fullfile(buoydat.buoylibdir, buoydat.SurgeFile);
            buoydat.ExcitationFile = fullfile(buoydat.buoylibdir, buoydat.ExcitationFile);
            buoydat.HydroCoeffsFile = fullfile(buoydat.buoylibdir, buoydat.HydroCoeffsFile);
        end
        
    end
    
    % set up the sea based on the spec
%     buoydat.SeaParameters = defaultseaparamaters(buoydat.SeaParameters, buoydat.BuoyParameters);

    % check the sea parameters are suitible for the buoy
    if max(buoydat.SeaParameters.sigma / (2*pi)) > buoydat.BuoyParameters.maxfreq
       
        warning('Max sea frequency greater than max buoy frequency, sea will be modified (higher frequencies stripped).')
        
        inds = find((buoydat.SeaParameters.sigma / (2*pi)) <= buoydat.BuoyParameters.maxfreq);
        
        buoydat.SeaParameters.sigma = buoydat.SeaParameters.sigma(inds);
        buoydat.SeaParameters.phase = buoydat.SeaParameters.phase(inds);
        buoydat.SeaParameters.amp = buoydat.SeaParameters.amp(inds);
        buoydat.SeaParameters.L = buoydat.SeaParameters.L(inds);
        buoydat.SeaParameters.wave_number = buoydat.SeaParameters.wave_number(inds);
        
        buoydat.SeaParameters.sigma_range(2) = max(buoydat.SeaParameters.sigma);
        
    end
    
    if min(buoydat.SeaParameters.sigma / (2*pi)) < buoydat.BuoyParameters.maxfreq
       
        warning('Min sea frequency less than min buoy frequency, sea will be modified (lower frequencies stripped).')
        
        inds = find((buoydat.SeaParameters.sigma / (2*pi)) >= buoydat.BuoyParameters.minfreq);
        
        buoydat.SeaParameters.sigma = buoydat.SeaParameters.sigma(inds);
        buoydat.SeaParameters.phase = buoydat.SeaParameters.phase(inds);
        buoydat.SeaParameters.amp = buoydat.SeaParameters.amp(inds);
        buoydat.SeaParameters.L = buoydat.SeaParameters.L(inds);
        buoydat.SeaParameters.wave_number = buoydat.SeaParameters.wave_number(inds);
        
        buoydat.SeaParameters.sigma_range(1) = min(buoydat.SeaParameters.sigma);
        
    end

    % Find the excitation forces
    [buoydat.BuoyParameters.heave_excit_force, buoydat.BuoyParameters.surge_excit_force] = ...
                unitampexcitforces(buoydat.SeaParameters.sigma, buoydat.ExcitationFile);
    
    % load up alphas and betas for this cylinder
    temp = load(buoydat.HeaveFile, 'alpha', 'beta');
    buoydat.BuoyParameters.Halpha = temp.alpha;
    buoydat.BuoyParameters.Hbeta = temp.beta;
    
    temp = load(buoydat.SurgeFile, 'alpha', 'beta');
    
    if isstruct(temp)
        buoydat.BuoyParameters.Salpha = temp.alpha;
        buoydat.BuoyParameters.Sbeta = temp.beta;
    else
        buoydat.BuoyParameters.Salpha = temp(:,1);
        buoydat.BuoyParameters.Sbeta = temp(:,2);
    end

%     % Find the number of alpha and betas used
%     num_coefficients = length(buoydat.BuoyParameters.Halpha);
%     num_coefficients_S = length(buoydat.BuoyParameters.Salpha);

    % Load the variables necessary to calculate the added mass for the
    % buoy, it is expected that this will be an ascii text data file
    hydrocoeffs = load(buoydat.HydroCoeffsFile, '-ascii');

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

    Hadded_mass = [(2 * pi ./ Hadded_mass(:,1)) , Hadded_mass(:,2) .* buoydat.BuoyParameters.rho];
    buoydat.BuoyParameters.HM_infinity = Hadded_mass(end);
    
    Sadded_mass = [(2 * pi ./ Sadded_mass(:,1)) , Sadded_mass(:,2) .* buoydat.BuoyParameters.rho];
    buoydat.BuoyParameters.SM_infinity = Sadded_mass(end);
        
    initial_vel_ext = 0.0;

    buoydat.BuoyParameters.velocity(:,1) = [initial_vel_ext, 0];
    
    % -----------------------------------------------------
    
    % set the number of radiation coefficients to use in the sim if it has
    % not already been specified
    if ~isfield(buoydat, 'NRadiationCoefs')
        buoydat.NRadiationCoefs = numel(buoydat.BuoyParameters.Halpha);
    else
        if buoydat.NRadiationCoefs < 1
            error('You must specify at least one radiation force coeficient.');
        elseif buoydat.NRadiationCoefs < 5
            warning('You have specified the use of less than 5 radiation coefficients, the accuracy may not be reliable');
        end
    end

    % construct the initial conditions of the system
    buoydat.ODESim.SolutionComponents.BuoyPositionHeave.InitialConditions = 0;
    buoydat.ODESim.SolutionComponents.BuoyPositionHeave.AbsTol = buoydat.BuoyParameters.draft / 20;
    
    buoydat.ODESim.SolutionComponents.BuoyVelocityHeave.InitialConditions = 0;
    buoydat.ODESim.SolutionComponents.BuoyVelocityHeave.AbsTol = 0.05;
    
    buoydat.ODESim.SolutionComponents.BuoyPositionSurge.InitialConditions = 0;
    buoydat.ODESim.SolutionComponents.BuoyPositionSurge.AbsTol = buoydat.BuoyParameters.a / 20;
    
    buoydat.ODESim.SolutionComponents.BuoyVelocitySurge.InitialConditions = 0;
    buoydat.ODESim.SolutionComponents.BuoyVelocitySurge.AbsTol = 0.05;
    
    % Handle the  the heave and surge radiation force component setup
    buoydat.ODESim.SolutionComponents.BuoyRadiationHeave.InitialConditions = zeros(1, buoydat.NRadiationCoefs);
    buoydat.ODESim.SolutionComponents.BuoyRadiationSurge.InitialConditions = zeros(1, buoydat.NRadiationCoefs);
    % Set the absolute tolerances on the hydrodynamic forces as the
    % force necessary to accelerate the mass of the buoy at a rate of
    % 0.01 m/s^2
    buoydat.ODESim.SolutionComponents.BuoyRadiationHeave.AbsTol = ...
        repmat(0.01 * buoydat.BuoyParameters.mass_external, 1, buoydat.NRadiationCoefs);
    buoydat.ODESim.SolutionComponents.BuoyRadiationSurge.AbsTol = ...
        repmat(0.01 * buoydat.BuoyParameters.mass_external, 1, buoydat.NRadiationCoefs);
        
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
    buoydat.EndStopks = (buoydat.BuoyParameters.rho * ...
                         buoydat.BuoyParameters.g * ...
                         pi * buoydat.BuoyParameters.a^2) / f;
                     
    
    % Set heave and surge not to be constrained if not set already
	buoydat.SeaParameters = setfieldifabsent(buoydat.SeaParameters, 'ConstrainSurge', 1);
    buoydat.SeaParameters = setfieldifabsent(buoydat.SeaParameters, 'ConstrainHeave', 1);
    
end