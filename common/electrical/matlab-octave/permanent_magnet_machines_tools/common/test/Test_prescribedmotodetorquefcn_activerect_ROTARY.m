% Test_prescribedmotodetorquefcn_activerect_ROTARY.m
%
% This script sets up a simulation of a radial flux permanent magnet
% machine with an average model of an active rectifier, using the full
% three-phase detailed representation of the machine, and controlled using
% field oriented control.
%

%% Set up the design

% Before we start, we'll get rid of anything in the workspace that might
% confuse us, i.e. any existing design and simoptions structure
clear design simoptions

design.Poles = 12;
% Choose the number of Phases, the conventional 3
design.Phases = 3;
% The desired number of layers is stored in 'CoilLayers'
design.CoilLayers = 2;
design.Qc = design.Phases * design.Poles;
design.qc = fr(design.Qc,design.Poles*design.Phases);
% Specify the actual coil slot pitch in slots
design.yd = 4;
design.CoilFillFactor = 0.6;
design.CoilTurns = 100;
design.Branches = 1;

% The outer radius of the rotor back iron
design.Ryo = 95e-3;
% The height of the magnets in the radial direction
design.tm = 6.3e-3;
% The thickness of the rotor back iron
design.tbi = 29.9523e-3;
% The thickness of the stator yoke (the part which the slots sit on)
design.ty = 17.4e-3;
% The radial height of a slot opening
design.tc = 16.4960e-003;
design.tc(2) = 2.1995e-003;
% The height of the slot shoe at the point it joins the slot wall
design.tsb = 2.604e-3;
% The height of the slot shoe at the point it joins the slot wall
design.tsg = 1.7364e-3;
% The radial air gap length between rotor and stator
design.g = 2e-3;
% magnet pitch in radians
design.thetam = (tau / design.Poles) * 0.8;
% coil slot opening pitch in radians ( space a coil has to fit in a slot )
% at the slot end closest to the air gap
design.thetacg = 0.5 * tau ()/design.Qc; % 84.7661e-003;
% coil slot opening pitch in radians ( space a coil has to fit in a slot )
% at the slot end closest to the armeture yoke
design.thetacy = design.thetacg;
% the pitch in radians of gap between the shoes that partially cover the
% slot openings
design.thetasg = 0.5 * design.thetacg;
% stack length of the machine, the depth along the cylindrical axis of the
% machine
design.ls = 88.9e-3;

design.MagnetSkew = 0.33;
design.NSkewMagnetsPerPole = 10;
design.ArmatureType = 'external';
design.MagnetPolarisation = 'radial';

% We must also specify some materials in the machine
design.MagFEASimMaterials.AirGap = 'Air';
design.MagFEASimMaterials.Magnet = 'NdFeB 40 MGOe';
design.MagFEASimMaterials.FieldBackIron = '1117 Steel';
design.MagFEASimMaterials.ArmatureYoke = design.MagFEASimMaterials.FieldBackIron;
design.MagFEASimMaterials.ArmatureCoil = '36 AWG';

% create the empty simoptions structure
simoptions.LoadModel = 'Machine Side Power Converter';
% complete the design
design = completedesign_RADIAL_SLOTTED(design, simoptions);

% converter parameters
RPM = 1800 / (design.Poles/2);
design.MachineSidePowerConverter.Rds = 0.2e-3;
% design.MachineSidePowerConverter.Vdc = 1.1 * (3 * sqrt(6) / pi) * ...
%     peakemfest_ROTARY (design.FluxLinkagePhasePeak, rpm2omega (RPM), design.Poles/2) / sqrt(2);
design.MachineSidePowerConverter.Vdc = 800;

%% Simulating the Machine

% GetVariableGapForce is a flag which determines whether
% simfun_RADIAL_SLOTTED does a series of force simulations at progressively
% smaller air gaps to create data on how the force varies with respect to
% this. We set this to false to save time.
simoptions.GetVariableGapForce = false;

% now run simfun_RADIAL_SLOTTED, this function performs various simulations
% of the machine
[design, simoptions] = simfun_RADIAL_SLOTTED(design, simoptions);

simoptions.Evaluation = designandevaloptions_RADIAL_SLOTTED();

% Now we are ready to post-process
[design, simoptions] = finfun_RADIAL_SLOTTED (design, simoptions);

%% Dynamic Simulation
%

TSpan = [0 10];
FOCStart = min ([5, TSpan(2) / 2]);

% Rectifier Set up

% following work for Vdc = 4000
% Kp_d = 2;
% Ki_d = 10;
% Kd_d = 0;
% Kp_q = 2;
% Ki_q = 10;
% Kd_q = 0;

% Kp_d = 0.731;
% Ki_d = 220;
% Kd_d = 0;
% Kp_q = Kp_d;
% Ki_q = Ki_d;
% Kd_q = Kd_d;

% set Vdc to be twice the natural rectification voltage
% 
% design.MachineSidePowerConverter.Rds = 0.2e-3;
% design.MachineSidePowerConverter.REquivalent = eye (3) .* (design.PhaseResistance + design.MachineSidePowerConverter.Rds);

% design.FOControl.PI_d = pidController ( Kp_d, Ki_d, Kd_d, 'MaxOut', design.MachineSidePowerConverter.Vdc, 'MinOut', -design.MachineSidePowerConverter.Vdc, 'InitialTime', FOCStart );
% design.FOControl.PI_q = pidController ( Kp_q, Ki_q, Kd_q, 'MaxOut', design.MachineSidePowerConverter.Vdc, 'MinOut', -design.MachineSidePowerConverter.Vdc, 'InitialTime', FOCStart );


simoptions.Isdref = 0; % flux reference (should always be zero for a PM machine)
simoptions.Isqref = -11.66; % torque reference (sets torque)
simoptions.FOCTApp = FOCStart;

simoptions = simsetup_ROTARY ( design, [], [], ...
                               'EvalFcn', 'prescribedmotodetorquefcn_activerect_ROTARY', ...
                               'Rpm', RPM, ...
                               'TSpan', TSpan, ...
                               'RampPoles', 20, ...
                               'simoptions', simoptions, ...
                               'MinPointsPerPole', 10 );

simoptions.ODESim.Solver = 'ode.rkfixed';

simoptions.ODESim.MaxStep = 1e-3;

% prescribedmotfinfun_ROTARY set up  
simoptions.ODESim.PostPreProcFcn = 'prescribedmotfinfun_RADIAL_SLOTTED';
simoptions.ODESim.MaxStep = min ([simoptions.ODESim.MaxStep, design.FOControl.MaxTimeStep]);

% With the simulation options set up in the simoptions structure we can use
% the common simulation function simulatemachine_AM to actually perform the
% simulation, and post-processes the results
simoptions.reltol = 1e-6;
simoptions.ODESim.OutputFcn = 'Test_prescribedmotodetorquefcn_activerect_ROTARY_ouputfcn';

[T, Y, results, design, simoptions] = simulatemachine_AM ( design, ...
                                                           simoptions );
                                                       
% The main output of the simulation is the phase currents
figure; 
yyaxis left
plot(T, results.EMF);
ylabel('Phase EMFs (V)');
yyaxis right
plot(T, Y);
ylabel('Phase Currents (A)');
xlabel('Time (s)');

figure; 
yyaxis left
plot(T, results.Vabc);
ylabel('Vabc (V)');
yyaxis right
plot(T, Y);
ylabel('Phase Currents (A)');
xlabel('Time (s)');

figure; 
plot (T, results.Tqpto);
xlabel('Time (s)');
ylabel('Torque (Nm)');


figure; 
% yyaxis left
plot (T, results.Idq);
hold on
plot (T([1,end]), [simoptions.Isqref; simoptions.Isqref], ':k');
xlabel('Time (s)');
ylabel('Current');
% yyaxis right
% plot (T, results.Vabc);
% ylabel('Voltage');
% hold off
% legend ('Id', 'Iq', 'I0', 'Isqref', 'Va', 'Vb', 'Vc');
legend ('Id', 'Iq', 'I0', 'Isqref');


figure; 
yyaxis left
plot (T, results.PI_vals(:,1:2));
xlabel('Time (s)');
ylabel('Error');
yyaxis right
plot (T, results.PI_vals(:,3:4));
ylabel('Integral');



