function simoptions = buoysimoptions(varargin)

    % the maximum allowed rms current density in a coil
    simoptions.maxAllowedJrms = 6e6;
    % the maximum allowed peak current density in the coil
    simoptions.maxAllowedJpeak = 10e6;
    % the minimum allowed rms voltage produced in a coil
    simoptions.minAllowedRMSEMF = 200;
    % the maximum allowed voltage produced by a coil
    simoptions.maxAllowedEMFpeak = [];

    % set the other simulation parameters 

    % The maximum allowed translator length
    simoptions.maxAllowedTLength = inf;
    % Maximum allowed translator displacement
    simoptions.maxAllowedxT = inf;
    % Maximum allowed armature displacement
    simoptions.maxAllowedxA = inf;
    % determines method used to calculate inductance
    simoptions.Lmode = 1;
    % the number of calculations to skip when producing output after the ode
    % solver finishes
    simoptions.skip = 1;
    % the time span of the simulation
    simoptions.ODESim.TimeSpan = [0, 60];   
    % Additional absolute tolerances on the components of the solution
    simoptions.ODESim.AbsTol = [];
    
end