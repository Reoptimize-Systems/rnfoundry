% mpga_machine_common_setup_options

mpgaoptions.GGAP = .8;           % Generation gap, how many new individuals are created
mpgaoptions.INSR = .9;           % Insertion rate, how many of the offspring are inserted
mpgaoptions.XOVR =  1;           % Crossover rate
mpgaoptions.SP = 2;              % Selective Pressure
mpgaoptions.MUTR = 1;            % Mutation rate; only a factor;
mpgaoptions.MIGR = 0.2;          % Migration rate between subpopulations
mpgaoptions.MIGGEN = 10;         % Number of generations between migration (isolation time)

mpgaoptions.TERMEXACT = 1e-4;    % Value for termination if minimum reached

mpgaoptions.SEL_F = 'sus';       % Name of selection function
mpgaoptions.XOV_F = 'recint';    % Name of recombination function for individuals
mpgaoptions.MUT_F = 'mutbga';    % Name of mutation function

% set up the penalties, empty penalties will be ignored

% the maximum allowed rms current density in a coil
simoptions.maxAllowedJrms = 6e6;
% the maximum allowed peak current density in the coil
simoptions.maxAllowedJpeak = 10e6;
% the minimum allowed rms voltage produced in a coil
simoptions.minAllowedRMSEMF = 200;
% the maximum allowed voltage produced by a coil
simoptions.maxAllowedEMFpeak = [];

% set the other simulation parameters 

% The maximum allowed translator length, this is a hard limit, not
% determined by a penalty. The number of Poles in the design will be
% modified if exceeded
simoptions.maxAllowedTLength = 5;
% The maximum allowed translator displacement. An end-stop will be applied
% if this limit is reached
simoptions.maxAllowedxT = inf;
% determines method used to calculate inductance
simoptions.Lmode = 1;
% the initial currents in the coils at t=0
simoptions.ODESim.InitialConditions = [0, 0, 0];
% the number of calculations to skip when producing output after the ode
% solver finishes
simoptions.ODESim.ResultsTSkip = 1;
% the time span of the simulation
simoptions.ODESim.TimeSpan = [0, 500];   
% Additional absolute tolerances on the components of the solution
simoptions.abstol = [];
% The peak frequency of a PM Specturm
SeaParameters.peak_freq = 1/9;
% The range of frequencies to be included in the spectrum
SeaParameters.sigma_range = [2*pi*0.055, 0.1*(0.6)^0.5 + 2.24];
% The number of frequencies to be produced
SeaParameters.freqcount = 50;
% The depth of the water
SeaParameters.water_depth = 50;
% Create the sea
simoptions.BuoySim.SeaParameters = defaultseaparamaters(SeaParameters);
% set the tether length
simoptions.BuoySim.tether_length = 6;
% Choose a size for the wave farm
simoptions.FarmSize = 10e6;
% Choose a buoy to use in the simulation
simoptions.buoynum = 37;

