function [T, Y, results, design, params] = ddlinearsystemsim(design, simoptions, simfunction, varargin)
% Finds all the forces for the cylindrical buoy and uses the matlab ode
% solvers to simulate the system together with the snapper machine model
%
% Syntax
%
% [T, Y, output, results, params] = ddlinearsystemsim(design, simoptions)
% 
% [T, Y, output, results, params] = ddlinearsystemsim(design, ...
%                                     simoptions, 'ParameterName', ...
%                                     ParameterValue, ...)
%     
% Description
%
% ddlinearsystemsim(design, simoptions, simfunction) simulates the combined
% bouy and direct-drive linear machine WEC system described by the machine
% design described in the 'design' structure and the simulation options in
% the simoptions structure.
%
% design is a structure containing information for the machine simulation:
%
%
% simoptions is another structure containing some options for the
% simulation. It can have the following fields:
%
%   tspan - (1 x 2) vector of time values between which to simulate the
%           machine.
%
%   IC - initial conditions for the snapper electrical simulation,
%   typically [0 0 0 0 0]
%
%   skip - optional value which determines how many of the internally
%   calculated values of the ode function are recalculated after the system
%   is solved. e.g if skip were 2, every other point in T output by the ode
%   would be recalculated, if 3, every third point and so on
%
% wholesystemsim_Snapper(..., 'ParameterName', ParameterValue, ...) also
% allows the specification of some additional parameters describing the
% buoy. ParameterName options are as follows:
%
%   'HeaveFile' - a string containing the file name of the data file
%   containing the coefficients to be used for the buoy, if not supplied,
%   the default 'coefficients_full_900_20.mat' is used
%
%   'SurgeFile' - a string containing the file name of the data file
%   containing the stuff!! to be used. If not supplied, the default
%   'L_kernal_d121009_full_cyl_surge.mat' is used.
%
%   'BuoyVarsFile' - a string containing the file name of the data file
%   containing the other stuff!! to be used. If not supplied, the default
%   'fullsize_double_cyl.1' is used. 
%
%   'BuoyParameters' - a structure containing the following fields
%   describing the buoy and sea conditions:
%
%
%     rho: the density of water, default 1025;
%
%     g: the acceleration due to gravity, default 9.81;
%
%     a: m radius of cylinder, default 4.6;
%
%     mass_external: Weight of buoypi * 4.6^2 * 10 * 30; 
%
%     water_depth = water depth, default 30;
%     
%     % 
%     sigma: Single Freq MONo sea, default 2*pi*0.10;
%
%     amp: m wave amplitude, default 1;
%     
%     phase: wave phasse, default pi/2;
%
%
% Examples
%
% 
% [T, Y, output, results, params] = ddlinearsystemsim(design, simoptions)
% 
% [T, Y, output, results, params] = ddlinearsystemsim(design, simoptions,...
%                                  'SurgeFile', 'L_kernal_d121009_full_cyl_surge_New.mat',...
%                                  'HeaveFile', 'coefficients_full_900_20_New.mat.mat')
%                              
% Output
%
% TODO: fill in output descriptions
                             


    % create an input parser to parse the optional arguments that can be
    % passed in, mostly these are things such as the file names of the
    % cylinder files. This allows us to enter the file names as optional
    % name-value pairs, using defaults if not present. e.g.
    %
    %   [...] = ddlinearsystemsim(design, simoptions, 'BuoyVarsFile', 'fullsize_double_cyl_New.1')
    %
    % loads the 'fullsize_double_cyl_New.1' file instead of the default. As
    % the HeaveFile and SurgeFile values are not set,
    % 'coefficients_full_900_20.mat' and
    % 'L_kernal_d121009_full_cyl_surge.mat' are loaded to get the alpha and
    % beta values
    
    ip = inputParser;
    ip.addRequired('design', @isstruct);
    ip.addRequired('simoptions', @isstruct);
    ip.addRequired('simfunction', @(x) isa(x, 'function_handle'));
    ip.addParamValue('HeaveFile', 'heave_coefficients_cyl_2di_1dr_d020610.mat', @ischar);
    ip.addParamValue('SurgeFile', 'surge_coefficients_cyl_2di_1dr.mat', @ischar);
    ip.addParamValue('BuoyVarsFile', 'cyl_d3103v4.1', @ischar);
    ip.addParamValue('HydroCoeffsFile', 'cyl_d3103v4.2', @ischar)
    ip.addParamValue('BuoyParameters', defaultbuoyparameters, @isstruct);
    ip.addParamValue('SeaParameters', defaultseaparamaters, @isstruct);
    ip.addParamValue('LegacyMode', false, @isscalar);
    ip.addParamValue('RootDir', '', @ischar);
    ip.FunctionName = 'POSITION_RICHARD';
    
    % Parse the input arguments, the optional arguments are then stored as
    % members of a structure in ip.Results
    ip.parse(design, simoptions, simfunction, varargin{:});
    
    % Get the correct file locations based on their relative paths to the
    % file wholesystemsim_Snapper.m. This allows us to easily farm out jobs
    % to Unix file systems
    HeaveFile = fullfile(ip.Results.RootDir, ip.Results.HeaveFile);
    SurgeFile = fullfile(ip.Results.RootDir, ip.Results.SurgeFile);
    BuoyVarsFile = fullfile(ip.Results.RootDir, ip.Results.BuoyVarsFile);
    HydroCoeffsFile = fullfile(ip.Results.RootDir, ip.Results.HydroCoeffsFile);
    
    % Get some initial buoy parameters, if not supplied, these are
    % generated by the 'defaultbuoyparamters' function
    params = ip.Results.BuoyParameters;
    
    seaparams = ip.Results.SeaParameters;

    % ************************************************************************

    % Find the excitation forces
    [params.excit_force, params.surge_force] = wamit_excitation(seaparams.sigma, HydroCoeffsFile);
    
    % load up alphas and betas for this cylinder
    temp = load(HeaveFile, 'alpha', 'beta'); % full size
    params.Halpha = temp.alpha;
    params.Hbeta = temp.beta;
    
    temp = load(SurgeFile, 'alpha', 'beta');
    if isstruct(temp)
        params.Salpha = temp.alpha;
        params.Sbeta = temp.beta;
    else
        params.Salpha = temp(:,1);
        params.Sbeta = temp(:,2);
    end

    % Find the number of alpha and betas used
    num_coefficients = length(params.Halpha);
    num_coefficients_S = length(params.Salpha);

    % Load the variables necessary to calculate the added mass for the
    % buoy, ARE THESE ALWAYS STORED AS ASCII DATA? CHECK
    buoyvars = load(BuoyVarsFile, '-ascii');

    Hadded_mass =[0 0]; % [sigma, added mass]
    Sadded_mass =[0 0]; % [sigma, added mass]
    
    % for n = 1:length(exp_sized_cyl)
    %     if exp_sized_cyl(n,2) ==3 && exp_sized_cyl(n,3) ==3
    %         added_mass (end+1, 1) =exp_sized_cyl(n,1); % wave number
    %         added_mass (end, 2) = exp_sized_cyl(n,4); % dimensionlised added mass
    %     end
    % end

    for n = 1:length(buoyvars)
        if buoyvars(n,2) ==3 && buoyvars(n,3) ==3
            Hadded_mass (end+1, 1) = buoyvars(n,1); % wave number
            Hadded_mass (end, 2) = buoyvars(n,4); % dimensionlised added mass
        end
        if buoyvars(n,2) ==1 && buoyvars(n,3) ==1
            Sadded_mass (end+1, 1) = buoyvars(n,1); % wave number
            Sadded_mass (end, 2) = buoyvars(n,4); % dimensionlised added mass
        end
    end
    Hadded_mass(1,:) = [];% Removes zero first entry
    Sadded_mass(1,:) = [];% Removes zero first entry

    Hadded_mass = [(2 * pi ./ Hadded_mass(:,1)) , Hadded_mass(:,2) .* params.rho];
    params.HM_infinity = Hadded_mass(end);
    Sadded_mass = [(2 * pi ./ Sadded_mass(:,1)) , Sadded_mass(:,2) .* params.rho];
    params.SM_infinity = Sadded_mass(end);
        
    initial_vel_ext = 0.0;

    params.velocity(:,1) = [initial_vel_ext, 0];
    
    % -----------------------------------------------------

    I_H = zeros(num_coefficients,1);
    I_S = zeros(num_coefficients_S,1);

    initial_conditions = [I_H; I_S]';
    
    if ~isfield(simoptions, 'tspan')
        simoptions.ODESim.TimeSpan = [0,1];
    end

    % Now perform the simulation by calling the appropriate function for
    % the machine
    [T, Y, results, design] = feval(simfunction, design, simoptions, params, seaparams, initial_conditions);
    
    % Calculate the wave heights
    results.wave_height = zeros(length(T),1);

    for k = 1:length(T)
        time_vector = T(k) * ones(1, length(seaparams.phase));
        results.wave_height(k) = sum(real(seaparams.amp .* exp(-i .* (seaparams.sigma .* time_vector - seaparams.phase))));
    end
    
end