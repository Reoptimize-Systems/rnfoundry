% Test_slotlessfemmprob_torus

design.taupm = 1;
design.taumm = 0.8;

design.g = 0.1; 
design.tc = 0.15;
design.ty = 0.15;
design.tauco = 0.333 * design.taupm;
design.tauci = 0.7 * design.tauco;

design.tm = 0.15;
design.tbi = [0.1, 0.2];

design.Rmo = 10;
design.Rmi = 9;

design.CoilTurns = 100;

% Matlib = parsematlib_mfemm(fullfile(fileparts(which('mfemm_parsematlib.m')), 'matlib.dat'));

% FemmProblem.Materials = Matlib([1, 47, 2]);

design.MagnetMaterial = 'NdFeB 32 MGOe';
design.BackIronMaterial = '1117 Steel';
design.CoilMaterial = '36 AWG';

FemmProblem = slotlessfemmprob_torus(design, 'NStages', 2, 'DrawCoils', true);


filename = 'test.fem';

writefemmfile(filename, FemmProblem)

openfemm;

opendocument(fullfile(pwd, filename))