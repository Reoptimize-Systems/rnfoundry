% Test_simfun_TORUS_SLOTTED

% FemmProblem = newproblem_mfemm('planar');

% %   Materials
% Matlib = parsematlib_mfemm(fullfile(fileparts(which('mfemm_parsematlib.m')), 'matlib.dat'));
% 
% FemmProblem.Materials = Matlib([1, 47, 2]);


design.taupm = 1;
design.taumm = 0.8;

design.g = 0.1; 
design.tc = 0.15;
design.tsb = 0.025;
design.tsg = 0;%0.5 * design.tsb;
design.ty = 0.15;

design.Phases = 3;
design.Poles = 28;
design.Qc = 7;
design.Qs = design.Phases * 2 * design.Qc;
design.yd = 1;
slotsperpole = design.Qs / design.Poles;

% design.tauco = 0.333 * design.taupm;

design.tm = 0.15;
design.tbi = [0.1, 0.2];

design.Rmo = 10;
design.Rmi = 9;

design.NPhaseCoils = design.Qc;

design.taucs = 0.8 * design.taupm / slotsperpole;
design.tausgm = design.taucs * 0.3;

design.Dc = design.taumm / 100;
design.CoilFillFactor = 0.7;

design.CoilLayers = 1;
design.Hc = design.tc / design.CoilLayers;
design.CoilTurns = 250;

design.MagnetMaterial = 'NdFeB 32 MGOe';
design.BackIronMaterial = '1117 Steel';
design.YokeMaterial = design.BackIronMaterial;
design.CoilMaterial = '36 AWG';

design.RgVRc = 10;

simoptions = struct();
simoptions.GetVariableGapForce = false;

[design, simoptions] = simfun_TORUS_SLOTTED(design, simoptions);

% openprobleminfemm_mfemm(design.FemmProblem);

[design, simoptions] = finfun_TORUS_SLOTTED(design, simoptions);

% plotslm(design.slm_psidot)


