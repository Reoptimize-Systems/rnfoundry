% +example_rm3/run_wecsim.m
%
% Example script demonstrating the simulation of the Reference Model 3
% buoy/spar/plate system with a basic PTO
%
%

% ensure some variables are cleared between runs (even though it's probably
% not that important)
clear waves simu hsys mbsys float_hbody spar_hbody hydro_mbnodes hydro_mbbodies hydro_mbelements lgsett wsobj

%% Hydro Simulation Data
simu = wsim.simSettings (fullfile ( wecsim_rootdir(), '..', '..', 'doc', 'sphinx', 'examples', '+example_rm3'));  % Create the Simulation Variable
simu.startTime = 0;                   % Simulation Start Time [s]
simu.endTime=400;                     % Simulation End Time [s]
simu.dt = 0.1;                        % Simulation time-step [s]
simu.rampT = 100;                     % Wave Ramp Time Length [s]
simu.multibodySolver = 'MBDyn';
simu.b2b = true;                      % BOdy-to-body interaction

%% Wave Information 

% noWaveCIC, no waves with radiation CIC  
% waves = waveClass('noWaveCIC');       %Create the Wave Variable and Specify Type      
%%%%%%%%%%%%%%%%%%%

% Regular Waves
waves = wsim.waveSettings ('regularCIC');        %Create the Wave Variable and Specify Type                               
waves.H = 2.5;                          %Wave Height [m]
waves.T = 8;                            %Wave Period [s]
%%%%%%%%%%%%%%%%%%

% Regular Waves
% waves = wsim.waveSettings ('regular');        %Create the Wave Variable and Specify Type                               
% waves.H = 2.5;                          %Wave Height [m]
% waves.T = 8;                            %Wave Period [s]
% simu.ssCalc = 1;
%%%%%%%%%%%%%%%%%%%

% Irregular Waves using PM Spectrum with Convolution Integral Calculation
% waves = wsim.waveSettings ('irregular');       %Create the Wave Variable and Specify Type
% waves.H = 2.5;                        %Significant Wave Height [m]
% waves.T = 8;                          %Peak Period [s]
% waves.spectrumType = 'PM';
%%%%%%%%%%%%%%%%%%%

% % Irregular Waves using BS Spectrum with Convolution Integral Calculation
% waves = wsim.waveSettings ('irregular');       %Create the Wave Variable and Specify Type
% waves.H = 2.5;                        %Significant Wave Height [m]
% waves.T = 8;                          %Peak Period [s]
% waves.spectrumType = 'BS';
% %%%%%%%%%%%%%%%%%%%

% Irregular Waves using BS Spectrum with State Space Calculation
% waves = wsim.waveSettings ('irregular');       %Create the Wave Variable and Specify Type
% waves.H = 2.5;                        %Significant Wave Height [m]
% waves.T = 8;                          %Peak Period [s]
% waves.spectrumType = 'BS';
% simu.ssCalc = 1;	
%Control option to use state space model 
%%%%%%%%%%%%%%%%%%%

% Irregular Waves using User-Defined Spectrum
% waves = wsim.waveSettings ('irregularImport');  %Create the Wave Variable and Specify Type
% waves.spectrumDataFile = 'ndbcBuoyData.txt';  %Name of User-Defined Spectrum File [2,:] = [omega, Sf]
%%%%%%%%%%%%%%%%%%%

% User-Defined Time-Series
% waves = wsim.waveSettings ('userDefined');     %Create the Wave Variable and Specify Type
% waves.etaDataFile = 'umpqua46229_6_2008.mat'; % Name of User-Defined Time-Series File [:,2] = [time, wave_elev]
%%%%%%%%%%%%%%%%%%%

%% Hydrodynamic body system
%
% Sets up the bodies which interact with the waves (and each other)
%
%

% Float
float_hbody = wsim.hydroBody('float.mat', 'CaseDirectory', simu.caseDir);      
    %Create the wsim.hydroBody(1) Variable, Set Location of Hydrodynamic Data File 
    %and Body Number Within this File.   
float_hbody.mass = 'equilibrium';                   
    %Body Mass. The 'equilibrium' Option Sets it to the Displaced Water 
    %Weight.
float_hbody.momOfInertia = [20907301, 21306090.66, 37085481.11];  %Moment of Inertia [kg*m^2]     
float_hbody.geometryFile = 'float.stl';    %Location of Geomtry File

% Spar/Plate
spar_hbody = wsim.hydroBody('spar.mat', 'CaseDirectory', simu.caseDir); 
spar_hbody.mass = 'equilibrium';                   
spar_hbody.momOfInertia = [94419614.57, 94407091.24, 28542224.82];
spar_hbody.geometryFile = 'plate.stl'; 

% octave can't make an array of objects with square brackets (like
% [float_hbody, spar_hbody]
obj_array(1) = float_hbody;
obj_array(2) = spar_hbody;
% make a hydrosys object for simulation
hsys = wsim.hydroSystem (waves, simu, obj_array);

% set up transient simulation
hsys.initialiseHydrobodies ();
hsys.odeSimSetup ();

% generate the nodes and elements for simulation of the hydrodynamic system
% in MBDyn. One node and one body element are created for each hydrodynamic
% body interacting with the waves. These are then used in the larger
% multibody system with other joints, constraints and bodies
[hydro_mbnodes, hydro_mbbodies, hydro_mbelements] = hsys.makeMBDynComponents ();

%% Multibody dynamics system specification (mbdyn)

problem_options.ResidualTol = 1e-5;
problem_options.MaxIterations = 200;
problem_options.Output = {}; % 'iterations', 'solution', 'jacobian matrix', 'matrix condition number', 'solver condition number'
% problem_options.NonLinearSolver = mbdyn.pre.newtonRaphsonSolver ();
problem_options.NonLinearSolver = [];
% problem_options.LinearSolver = mbdyn.pre.linearSolver ('naive');
problem_options.LinearSolver = [];

[mbsys, initptodpos] = example_rm3.make_multibody_system ( waves, ...
                                                           simu, ...
                                                           hydro_mbnodes, ...
                                                           hydro_mbbodies, ...
                                                           hydro_mbelements, ...
                                                           problem_options );
                     
% draw it in a figure
mbsys.draw ( 'Mode', 'solid', ...
             'Light', true, ...
             'AxLims', [-30, 30; -30, 30; -35, 35], ...
             'Joints', false, ...
             'StructuralNodes', false)

% create the path where the MBDyn input file will be generated
mbdpath = fullfile (simu.caseDir, 'RM3.mbd');

%% Set up Power Take-Off (PTO)

% set up a simple linear spring-damper power take-off force based on a
% matlab function. See help for
k = 0;
c = 1200000;

forcefcn = @(time, xRpto, vRpto) -k*xRpto -c*vRpto;

% create a power take-off object attached to the two hydro nodes, with the
% force being based on the relative velocity and displacement along axis 3
% of the first (reference) node (in it's coordinate frame).
pto = wsim.linearPowerTakeOff ( hydro_mbnodes{2}, hydro_mbnodes{1}, 3, forcefcn );

%% Run the simulation

lgsett = wsim.loggingSettings ();

lgsett.positions = true;
lgsett.velocities = true;
lgsett.accelerations = true;
lgsett.nodeForcesAndMoments = true;
lgsett.nodeForcesAndMomentsUncorrected = true;
lgsett.forceHydro = true;
lgsett.forceExcitation = true;
lgsett.forceExcitationRamp = true;
lgsett.forceExcitationLin = true;
lgsett.forceExcitationNonLin = true;
lgsett.forceRadiationDamping = true;
lgsett.forceRestoring = true;
lgsett.forceMorrison = true;
lgsett.forceViscousDamping = true;
% lssett.ForceAddedMassUncorrected = false;
lgsett.forceAddedMass = true;
        
% create the wesim object
wsobj = wsim.wecSim ( hsys, mbsys, ...
                      'PTO', pto, ... % multiple PTOs may be added
                      'LoggingSettings', lgsett );

% initialise the simulation
wsobj.prepare ();

% run it and get the output data
datalog = wsobj.run ('TimeExecution', true);

%% Plot some results

datalog.plotVar ('Positions');

datalog.plotVar ('Velocities');
