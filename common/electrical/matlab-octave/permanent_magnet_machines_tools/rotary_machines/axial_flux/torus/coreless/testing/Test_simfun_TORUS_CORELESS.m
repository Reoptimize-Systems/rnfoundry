% Test_simfun_TORUS_CORELESS

design.g = 0.1; 

design.tc = 0.3;
    
design.taupm = 1;

design.Phases = 3;
design.taupcg = (4/3) * design.taupm * design.Phases;
design.tauco = 0.32 * design.taupcg;
design.tauci = 0.7 * design.tauco;
design.CoilTurns = 250;

design.taumm = 0.8;

design.tm = 0.15;

design.tbi = 0.05;

design.Rmo = 10;
design.Rmi = 9;
% Inner radius of coil
design.Rci = design.Rmi;
% Outer Radius of coil
design.Rco = design.Rmo;

design.Hc = design.tc;
design.Wc = (design.tauco - design.tauci) / 2;

design.Dc = design.taumm / 100;
design.CoilFillFactor = 0.7;

% Matlib = parsematlib_mfemm(fullfile(fileparts(which('mfemm_parsematlib.m')), 'matlib.dat'));

% FemmProblem.Materials = Matlib([1, 47, 2]);

design.MagFEASimMaterials.Magnet = 'NdFeB 32 MGOe';
design.MagFEASimMaterials.FieldBackIron = '1117 Steel';
design.MagFEASimMaterials.ArmatureCoil = '36 AWG';

design.WindingType = 'nonoverlapping';

simoptions = struct;

[design, simoptions] = simfun_TORUS_CORELESS(design, simoptions);






