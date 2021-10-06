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

Rphase = 0.7250;
Lphase = 0.0091;
FluxPhasePeak = 2.6281;
Poles = 32;
Qc = 24;

[design, simoptions] = machineprops_to_design_ROTARY ( Rphase, ...
                                                       Lphase, ...
                                                       Poles, ...
                                                       Qc, ...
                                                       'FluxPhasePeak', FluxPhasePeak, ...
                                                       'yd', 1, ...
                                                       'CoilLayers', 2 );

% create the empty simoptions structure
simoptions.LoadModel = 'Machine Side Power Converter';
% converter parameters
RPM = 100 / (design.Poles/2);
design.MachineSidePowerConverter.Rds = 0.2e-3;
% design.MachineSidePowerConverter.Vdc = 1.1 * (3 * sqrt(6) / pi) * ...
%     peakemfest_ROTARY (design.FluxLinkagePhasePeak, rpm2omega (RPM), design.Poles/2) / sqrt(2);
design.MachineSidePowerConverter.Vdc = 675;

% design.FOControl.DirectCurentKp = 32.7483; 
% design.FOControl.DirectCurentKi = 3.9301e+04;
% design.FOControl.QuadratureCurentKp = design.FOControl.DirectCurentKp;
% design.FOControl.QuadratureCurentKi = design.FOControl.DirectCurentKi;

% Now we are ready to post-process
[design, simoptions] = finfun_ROTARY (design, simoptions);

%% Dynamic Simulation
%

TSpan = [0, 10];
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


simoptions.Isdref = -9; % flux reference (should always be zero for a PM machine)
simoptions.Isqref = 0; % torque reference (sets torque)
simoptions.FOCTApp = FOCStart;

simoptions = simsetup_ROTARY ( design, [], [], ...
                               'EvalFcn', 'prescribedmotode_activerect_ROTARY', ...
                               'Rpm', RPM, ...
                               'TSpan', TSpan, ...
                               'RampPoles', 20, ...
                               'simoptions', simoptions, ...
                               'MinPointsPerPole', 10, ...
                               'ForceAddPhaseCurrentODESolutionComps', true);

simoptions.ODESim.Solver = 'ode.rkfixed';

simoptions.ODESim.MaxStep = 5e-5;

% prescribedmotfinfun_ROTARY set up  
simoptions.ODESim.PostPreProcFcn = 'prescribedmotfinfun_RADIAL_SLOTTED';
simoptions.ODESim.MaxStep = min ([simoptions.ODESim.MaxStep, design.FOControl.MaxTimeStep]);

% With the simulation options set up in the simoptions structure we can use
% the common simulation function simulatemachine_AM to actually perform the
% simulation, and post-processes the results
simoptions.reltol = 1e-6;
simoptions.ODESim.OutputFcn = 'Test_prescribedmotodetorquefcn_activerect_ROTARY_ouputfcn';

tic
[T, Y, results, design, simoptions] = simulatemachine_AM ( design, ...
                                                           simoptions );
toc

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



