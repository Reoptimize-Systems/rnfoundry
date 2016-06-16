% Test_evaluatedesign_TORUS_SLOTLESS
clear design simoptions
design.Phases = 3;

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

design.MagFEASimMaterials.Magnet = 'NdFeB 32 MGOe';
design.MagFEASimMaterials.FieldBackIron = '1117 Steel';
design.MagFEASimMaterials.ArmatureCoil = '36 AWG';

simoptions = struct;
simoptions.GetVariableGapForce = false;


% setup simulation options
simoptions = simsetup_ROTARY(design, 'simfun_TORUS_SLOTLESS', 'finfun_TORUS_SLOTLESS', ...
                                'Velocity', 1, 'TSpan', [0,10]);
                            
simoptions.reltol = 1e-4;
simoptions.abstol = repmat(0.001, 1, design.Phases);
simoptions.maxstep = (simoptions.ODESim.TimeSpan(2) - simoptions.ODESim.TimeSpan(1)) / 10000;

[score, design, simoptions, T, Y, results] = ...
    evaluatedesign_TORUS_SLOTLESS(design, simoptions);

