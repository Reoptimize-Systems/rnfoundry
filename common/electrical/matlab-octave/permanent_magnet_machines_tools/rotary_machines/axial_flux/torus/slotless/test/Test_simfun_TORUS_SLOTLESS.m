% Test_simfun_TORUS_SLOTLESS
clear design
design.Phases = 3;
design.Poles = 24;
design.Qc = 7;
design.Qs = design.Phases * 2 * design.Qc;
design.yd = 1;
design.NPhaseCoils = design.Qc;

design.taupm = 1;
design.taumm = 0.8;

design.g = 0.1; 
design.tc = 0.15;
design.ty = 0.15;
design.taucsm = design.taupm / design.Phases;
design.taupcg = design.Phases * design.taucsm;
design.tauco = 0.95 * design.taucsm;

design.tm = 0.15;
design.tbi = [0.1, 0.2];

design.Rmo = 10;
design.Rmi = 9;

design.Hc = design.tc;
design.Wc = design.tauco;

design.Dc = design.taumm / 100;
design.CoilFillFactor = 0.7;

design.CoilTurns = 250;

% Matlib = parsematlib_mfemm(fullfile(fileparts(which('mfemm_parsematlib.m')), 'matlib.dat'));

% FemmProblem.Materials = Matlib([1, 47, 2]);

design.MagnetMaterial = 'NdFeB 32 MGOe';
design.BackIronMaterial = '1117 Steel';
design.CoilMaterial = '36 AWG';

design.RlVRp = 10;

simoptions = struct;
simoptions.GetVariableGapForce = false;

[design, simoptions] = simfun_TORUS_SLOTLESS(design, simoptions);

[design, simoptions] = finfun_TORUS_SLOTLESS(design, simoptions);

