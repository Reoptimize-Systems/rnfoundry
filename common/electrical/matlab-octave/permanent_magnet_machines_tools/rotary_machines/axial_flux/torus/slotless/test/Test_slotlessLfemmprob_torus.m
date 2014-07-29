% Test_slotlessLfemmprob_torus

design.taupm = 1;
design.taumm = 0.8;

design.g = 0.1; 
design.tc = 0.15;
design.ty = 0.15;
design.Phases = 3;
design.taucsm = design.taupm/design.Phases; 
design.tauco = 0.95 * design.taucsm;
design.taupcg = design.Phases * design.taucsm;
design.tm = 0.15;
design.tbi = [0.1, 0.2];

design.Rmo = 10;
design.Rmi = 9;

design.CoilTurns = 100;

% Matlib = parsematlib_mfemm(fullfile(fileparts(which('mfemm_parsematlib.m')), 'matlib.dat'));

% FemmProblem.Materials = Matlib([1, 47, 2]);

design.MagSimMaterials.Magnet = 'NdFeB 32 MGOe';
design.MagSimMaterials.FieldIron = '1117 Steel';
design.MagSimMaterials.CoilWinding = '36 AWG';

FemmProblem = slotlessLfemmprob_torus(design, 'NStages', 2);

filename = 'test.fem';

writefemmfile(filename, FemmProblem)

openfemm;

opendocument(fullfile(pwd, filename))