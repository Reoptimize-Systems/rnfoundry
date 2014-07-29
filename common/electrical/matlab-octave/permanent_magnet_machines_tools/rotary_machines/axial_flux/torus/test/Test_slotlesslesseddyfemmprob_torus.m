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

design.Poles = pi * 2 * mean([design.Rmo, design.Rmi]) / design.taupm;

design.CoilTurns = 100;

% Matlib = parsematlib_mfemm(fullfile(fileparts(which('mfemm_parsematlib.m')), 'matlib.dat'));

% FemmProblem.Materials = Matlib([1, 47, 2]);

design.MagSimMaterials.Magnet = 'NdFeB 32 MGOe';
design.HcMag = 915e3;
design.MagSimMaterials.FieldIron = '1117 Steel';
design.MagSimMaterials.CoilWinding = '36 AWG';

RPM = 1000;

FemmProblem = slotlesseddyfemmprob_torus(design, RPM, 'NStages', 1);


filename = 'test.fem';

writefemmfile(filename, FemmProblem)

openfemm;

opendocument(fullfile(pwd, filename))