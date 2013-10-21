% Test_rotgen_three_phase_PM_machine_mdl.m
%
%

clear data

data.design.StatorType = 'so';
data.design.Poles = 2;
data.design.Phases = 3;
data.design.CoilLayers = 2;
data.design.Qc = data.design.Phases * data.design.Poles;
if data.design.CoilLayers == 1
    data.design.Qs = data.design.Phases * 2 * data.design.Qc;
elseif data.design.CoilLayers == 2
    data.design.Qs = data.design.Phases * 1 * data.design.Qc;
end
data.design.yd = 4;
data.design.thetap = 2*pi/data.design.Poles;
data.design.thetam = data.design.thetap * 0.8;
data.design.thetac = (2*pi / data.design.Qs) * 0.85;
data.design.thetasg = data.design.thetac * 0.95;
data.design.tm = 0.01;
data.design.tbi = 0.01;
data.design.ty = 0.01;
data.design.tc = 0.03;
data.design.tsb = 0.01;
data.design.tsg = 0; %0.01;
data.design.g = 3/1000;
data.design.Rmo = 0.5;
data.design.Rmi = 0.5;
data.design.ls = 0.3;

if strcmp(data.design.StatorType, 'si')
    data.design.Rmo = 0.5;
    data.design.Rmi = data.design.Rmi - data.design.tm;
    data.design.Rmm = mean([data.design.Rmi, data.design.Rmo]);
    data.design.Rci = data.design.Rmo + data.design.g + data.design.tsb;
    data.design.Rco = data.design.Rci + data.design.tc;
    data.design.Rcm = mean([data.design.Rci, data.design.Rco]);
    data.design.Rbo = data.design.Rmi;
    data.design.Rbi = data.design.Rbo - data.design.tbi;
    data.design.Rbm = mean([data.design.Rbo, data.design.Rbi]);
    data.design.Ryi = data.design.Rco;
    data.design.Ryo = data.design.Rco + data.design.ty;
    data.design.Rym = mean([data.design.Ryi, data.design.Ryo]);
elseif strcmp(data.design.StatorType, 'so')
    data.design.Rmi = 0.5;
    data.design.Rmo = data.design.Rmi + data.design.tm;
    data.design.Rmm = mean([data.design.Rmi, data.design.Rmo]);
    data.design.Rco = data.design.Rmi - data.design.g - data.design.tsb;
    data.design.Rci = data.design.Rco - data.design.tc;
    data.design.Rcm = mean([data.design.Rci, data.design.Rco]);
    data.design.Rbi = data.design.Rmo;
    data.design.Rbo = data.design.Rbi + data.design.tbi;
    data.design.Rbm = mean([data.design.Rbo, data.design.Rbi]);
    data.design.Ryo = data.design.Rci;
    data.design.Ryi = data.design.Ryo - data.design.ty;
    data.design.Rym = mean([data.design.Ryi, data.design.Ryo]);
end

data.design.Dc = data.design.Rcm * data.design.thetac / 100;
data.design.CoilFillFactor = 0.7;

data.design.Hc = data.design.tc / data.design.CoilLayers;
data.design.CoilTurns = 250;

data.design.NCoilsPerPhase = data.design.Qc / data.design.Phases;

data.design.MagnetMaterial = 'NdFeB 32 MGOe';
data.design.BackIronMaterial = '1117 Steel';
data.design.YokeMaterial = data.design.BackIronMaterial;
data.design.CoilMaterial = '36 AWG';

data.design.RlVRp = 10;

data.simoptions = struct();
data.simoptions.GetVariableGapForce = false;

[data.design, data.simoptions] = simfun_RADIAL_SLOTTED(data.design, data.simoptions);

[data.design, data.simoptions] = finfun_RADIAL_SLOTTED(data.design, data.simoptions);

%% Now set up the simulink simulation

open_system('Test_rotgen_three_phase_pm_machine');

set_param('Test_rotgen_three_phase_pm_machine/Three Phase Rotary PM Machine/Three Phase Rotary PM Machine S-Function', 'UserData', data);

